#!/usr/bin/env python

import fastaIO
import sys
import os
import os.path
import subprocess as subp
from fastaIO import Vividict
from collections import OrderedDict
from operator import itemgetter

args = sys.argv[1:]

def usage():
    print """
    Usage: 
    parse_mummer_overlap_for_mix.py <mummer_overlap_tab_file> <assembly_fasta_file> <run_name> <scaffolding> <contig_in_file> <contig_self_file>

    This script parses a mummer/nucmer overlap output file in tabular format, finding the contigs wholely contained within other contigs. These are then removed from the assembly. If rescaffolding the contigs from the scaffolding in the original assesmblies is wanted, put 'yes', 'Yes', 'y', 'Y', or any word starting with "Y". Other values will not trigger scaffolding. If merging the outputs of previous split runs, combine the contig_in and contig_self files from all previous runs and provide their paths. 
    
    
    """
    sys.exit(-1)

if (len(args) != 4 and len(args) != 5 and len(args) != 6) or sys.argv[1] == '-h' or sys.argv[1] == '-help' or sys.argv[1] == '-H' or sys.argv[1] == '-Help' or sys.argv[1] == '--h' or sys.argv[1] == '--help':
    usage()

'''
Internal data formats

seen dictionary:
seen[ref_name][end of overlap relative to ref] = [len query extension, query_name, seq of extension, end of overlap relative to query, orientation of query to actual ref, actual ref_name to add extension to, which (ref or query) are combined, actual query_name, ori to actual query]

good dictionary:
good[best_query_contig] = [ref_contig, end of ref_contig extended, end of query_contig used in extension, orientation of query to final orientation (either just to ref if it's a single contig or to the orientation of merged contigs)]

'''

def process_single(seen, good, assembly, bad, report_list, contigs, contig_self, ref_start, ref_end, query_start, query_end, ref_align_len, query_align_len, percent_id, ref_len, query_len, ref_coverage, query_coverage, frame, strand, ref_name, query_name, tag):
    if ref_name not in seen:
        seen[ref_name]["start"] = []
        seen[ref_name]["end"] = []
        
    if query_coverage == "100.00" and int(query_align_len) == int(query_len):
        bad[query_name] = 1
        if "." in query_name:
            super_num, contig_num = query_name.split(".")
            contig_self[super_num][contig_num][ref_name] = strand
        assembly.pop(query_name, None)
        contigs.pop(query_name, None)
        report_list.append(ref_name + "\t" + query_name + "\t" + "covered_100")
            
    else:
        ref_end_dif = int(ref_len) - int(ref_end)
        ref_start_dif = int(ref_start)
        
        '''
    seen[ref_name][end of overlap relative to ref] = [len query extension, query_name, seq of extension, end of overlap relative to query, orientation of query to actual ref, actual ref_name to add extension to, which (ref or query) are combined, actual query_name, ori to actual query]
        '''
        
        if strand == "1" and (int(ref_start) < 24 or int(ref_end) > int(ref_len) - 23):
            end_dif = int(query_len) - int(query_end)
            start_dif = int(query_start)
            
            if end_dif > ref_end_dif:
                extra_end_seq = assembly[query_name][-(end_dif-ref_end_dif+1):]
                seen[ref_name]["end"].append([len(extra_end_seq), query_name, extra_end_seq, "start", strand, ref_name, percent_id, query_name, strand])
                report_list.append("Grabbing sequence in process_single. 1st " + query_name + " " + ref_name)
            if start_dif > ref_start_dif:
                extra_start_seq = assembly[query_name][:(start_dif-ref_start_dif)]
                seen[ref_name]["start"].append([len(extra_start_seq), query_name, extra_start_seq, "end", strand, ref_name, percent_id, query_name, strand])
                report_list.append("Grabbing sequence in process_single. 2nd"+ " " + query_name + " " + ref_name)
        elif strand == "-1" and (int(ref_start) < 24 or int(ref_end) > int(ref_len) - 23):
            end_dif = int(query_end)
            start_dif = int(query_len) - int(query_start)
            
            if end_dif > ref_end_dif:
                extra_end_seq = fastaIO.reverse_complement(assembly[query_name][:(end_dif - ref_end_dif)])
                seen[ref_name]["end"].append([len(extra_end_seq), query_name, extra_end_seq, "end", strand, ref_name, percent_id, query_name, strand])
                report_list.append("Grabbing sequence in process_single. 3rd" + " " + query_name + " " + ref_name)
            if start_dif > ref_start_dif:
                extra_start_seq = fastaIO.reverse_complement(assembly[query_name][-(start_dif - ref_start_dif + 1):])
                seen[ref_name]["start"].append([len(extra_start_seq), query_name, extra_start_seq, "start", strand, ref_name, percent_id, query_name, strand]) 
                report_list.append("Grabbing sequence in process_single. 4th" + " " + query_name + " " + ref_name)
        elif float(query_coverage) >= 98.00 and float(percent_id) >= 98.00:
            bad[query_name] = 1
            if "." in query_name:
                super_num, contig_num = query_name.split(".")
                contig_self[super_num][contig_num][ref_name] = strand
            assembly.pop(query_name, None)
            contigs.pop(query_name, None)
            report_list.append(ref_name + "\t" + query_name + "\t" + "covered_98")
    return bad, seen, assembly, good, contigs, contig_self

def process_ref_combined(seen, good, assembly, bad, report_list, contigs, contig_self, ref_start, ref_end, query_start, query_end, ref_align_len, query_align_len, percent_id, ref_len, query_len, ref_coverage, query_coverage, frame, strand, ref_name, query_name, tag):
        
    if ref_name not in seen:
        seen[ref_name]["start"] = []
        seen[ref_name]["end"] = []

    combined_ref_name = good[ref_name][0]
    combined_overlap_end = good[ref_name][1]
    ref_overlap_end = good[ref_name][2]
    ref_to_combined_strand = good[ref_name][3]
    
    if query_coverage == "100.00" and int(query_align_len) == int(query_len) and float(percent_id) >= 97.9:
        bad[query_name] = 1
        if "." in query_name:
            super_num, contig_num = query_name.split(".")
            if strand == "1":
                if ref_to_combined_strand == "1":
                    contig_self[super_num][contig_num][combined_ref_name] = "1"
                else:
                    contig_self[super_num][contig_num][combined_ref_name] = "-1"
            else:
                if ref_to_combined_strand == "1":
                    contig_self[super_num][contig_num][combined_ref_name] = "-1"
                else:
                    contig_self[super_num][contig_num][combined_ref_name] = "1"
        assembly.pop(query_name, None)
        contigs.pop(query_name, None)
        report_list.append(ref_name + "\t" + query_name + "\t" + "covered_100")
        return bad, seen, assembly, good
            
    
    ref_end_dif = int(ref_len) - int(ref_end)
    ref_start_dif = int(ref_start)
    if strand == "1" and (int(ref_start) < 24 or int(ref_end) > int(ref_len) - 23):
        end_dif = int(query_len) - int(query_end)
        start_dif = int(query_start)
        '''
        seen[ref_name][end of overlap relative to ref] = [len query extension, query_name, seq of extension, end of overlap relative to query, orientation of query to actual ref, actual ref_name to add extension to, which (ref or query) are combined, actual query_name, ori to actual query]
        '''
        
        if end_dif > ref_end_dif:
            if ref_overlap_end == "start" and ref_to_combined_strand == "1":
                extra_end_seq = assembly[query_name][(end_dif - ref_end_dif):]
                seen[ref_name]["end"].append([len(extra_end_seq), query_name, extra_end_seq, "start", "1", combined_ref_name, percent_id, query_name, strand])
                report_list.append("Grabbing sequence in process_ref_combined. 1st" + " " + query_name + " " + ref_name + " " + combined_ref_name)
            elif ref_overlap_end == "start" and ref_to_combined_strand == "-1":
                extra_end_seq = fastaIO.reverse_complement(assembly[query_name][(end_dif - ref_end_dif):])
                seen[ref_name]["start"].append([len(extra_end_seq), query_name, extra_end_seq, "start", "-1", combined_ref_name, percent_id, query_name, strand])
                report_list.append("Grabbing sequence in process_ref_combined. 2nd" + " " + query_name + " " + ref_name + " " + combined_ref_name)
                
        if start_dif > ref_start_dif:
            if ref_overlap_end == "end" and ref_to_combined_strand == "1":
                extra_start_seq = assembly[query_name][:-(start_dif - ref_start_dif)]
                seen[ref_name]["start"].append([len(extra_start_seq), query_name, extra_start_seq, "end", "1", combined_ref_name, percent_id, query_name, strand])
                report_list.append("Grabbing sequence in process_ref_combined. 3rd" + " " + query_name + " " + ref_name + " " + combined_ref_name)
            elif ref_overlap_end == "end" and ref_to_combined_strand == "-1":
                extra_start_seq = fastaIO.reverse_complement(assembly[query_name][:-(start_dif - ref_start_dif)])
                seen[ref_name]["end"].append([len(extra_start_seq), query_name, extra_start_seq, "start", "-1", combined_ref_name, percent_id, query_name, strand])
                report_list.append("Grabbing sequence in process_ref_combined. 4th" + " " + query_name + " " + ref_name + " " + combined_ref_name)
                
    elif strand == "-1" and (int(ref_start) < 24 or int(ref_end) > int(ref_len) - 23):
        end_dif = int(query_end)
        start_dif = int(query_len) - int(query_start)
        
        if end_dif > ref_end_dif:
            if ref_overlap_end == "start" and ref_to_combined_strand == "1":
                extra_end_seq = fastaIO.reverse_complement(assembly[query_name][:(end_dif - ref_end_dif)])
                seen[ref_name]["end"].append([len(extra_end_seq), query_name, extra_end_seq, "end", "-1", combined_ref_name, percent_id, query_name, strand])
                report_list.append("Grabbing sequence in process_ref_combined. 5th" + " " + query_name + " " + ref_name + " " + combined_ref_name)
            elif ref_overlap_end == "start" and ref_to_combined_strand == "-1":
                extra_end_seq = assembly[query_name][:(end_dif - ref_end_dif)]
                seen[ref_name]["start"].append([len(extra_end_seq), query_name, extra_end_seq, "end", "1", combined_ref_name, percent_id, query_name, strand])
                report_list.append("Grabbing sequence in process_ref_combined. 6th" + " " + query_name + " " + ref_name + " " + combined_ref_name)
                    
        if start_dif > ref_start_dif:
            if ref_overlap_end == "end" and ref_to_combined_strand == "1":
                extra_start_seq = fastaIO.reverse_complement(assembly[query_name][-(start_dif - ref_start_dif + 1):])
                seen[ref_name]["start"].append([len(extra_start_seq), query_name, extra_start_seq, "end", "-1", combined_ref_name, percent_id, query_name, strand])
                report_list.append("Grabbing sequence in process_ref_combined. 7th" + " " + query_name + " " + ref_name + " " + combined_ref_name)
            elif ref_overlap_end == "end" and ref_to_combined_strand == "-1":
                extra_start_seq = assembly[query_name][-(start_dif - ref_start_dif + 1):]
                seen[ref_name]["end"].append([len(extra_start_seq), query_name, extra_start_seq, "start", "1", combined_ref_name, percent_id, query_name, strand])
                report_list.append("Grabbing sequence in process_ref_combined. 8th" + " " + query_name + " " + ref_name + " " + combined_ref_name)
    elif float(query_coverage) >= 98.00 and float(percent_id) >= 98.00:
        bad[query_name] = 1
        if "." in query_name:
            super_num, contig_num = query_name.split(".")
            if strand == "1":
                if ref_to_combined_strand == "1":
                    contig_self[super_num][contig_num][combined_ref_name] = "1"
                else:
                    contig_self[super_num][contig_num][combined_ref_name] = "-1"
            else:
                if ref_to_combined_strand == "1":
                    contig_self[super_num][contig_num][combined_ref_name] = "-1"
                else:
                    contig_self[super_num][contig_num][combined_ref_name] = "1"
        assembly.pop(query_name, None)
        contigs.pop(query_name, None)
        report_list.append(ref_name + "\t" + query_name + "\t" + "covered_98")
    return bad, seen, assembly, good, contig_self

def process_query_combined(seen, good, assembly, bad, report_list, contigs, ref_start, ref_end, query_start, query_end, ref_align_len, query_align_len, percent_id, ref_len, query_len, ref_coverage, query_coverage, frame, strand, ref_name, query_name, tag):
    if ref_name not in seen:
        seen[ref_name]["start"] = []
        seen[ref_name]["end"] = []
    
    combined_query_name = good[query_name][0]
    combined_overlap_end = good[query_name][1]
    query_overlap_end = good[query_name][2]
    query_to_combined_strand = good[query_name][3]
            
    
    ref_end_dif = int(ref_len) - int(ref_end)
    ref_start_dif = int(ref_start)
    if strand == "1" and (int(ref_start) < 24 or int(ref_end) > int(ref_len) - 23):
        end_dif = int(query_len) - int(query_end)
        start_dif = int(query_start)
        
        '''
        seen[ref_name][end of overlap relative to ref] = [len query extension, query_name, seq of extension, end of overlap relative to query, orientation of query to actual ref, actual ref_name to add extension to, which (ref or query) are combined, actual query_name, ori to actual query]
        '''
        
        if end_dif > ref_end_dif:
            if query_overlap_end == "end" and query_to_combined_strand == "1":
                extra_end_seq = assembly[combined_query_name][(end_dif - ref_end_dif):]
                seen[ref_name]["end"].append([len(extra_end_seq), query_name, extra_end_seq, "start", "1", ref_name, percent_id, combined_query_name, query_to_combined_strand])
                report_list.append("Grabbing sequence in process_query_combined. 1st" + " " + query_name + " " + ref_name + " " + combined_query_name)
            elif query_overlap_end == "end" and query_to_combined_strand == "-1":
                extra_end_seq = fastaIO.reverse_complement(assembly[combined_query_name][:-(end_dif - ref_end_dif)])
                seen[ref_name]["end"].append([len(extra_end_seq), query_name, extra_end_seq, "start", "1", ref_name, percent_id, combined_query_name, query_to_combined_strand])
                report_list.append("Grabbing sequence in process_query_combined. 2nd" + " " + query_name + " " + ref_name + " " + combined_query_name)
        if start_dif > ref_start_dif:
            if query_overlap_end == "start" and query_to_combined_strand == "1":
                extra_start_seq = assembly[combined_query_name][:-int(query_align_len)]
                seen[ref_name]["start"].append([len(extra_start_seq), query_name, extra_start_seq, "end", "1", ref_name, percent_id, combined_query_name, query_to_combined_strand])
                report_list.append("Grabbing sequence in process_query_combined. 3rd" + " " + query_name + " " + ref_name + " " + combined_query_name)
            elif query_overlap_end == "start" and query_to_combined_strand == "-1":
                extra_start_seq = fastaIO.reverse_complement(assembly[combined_query_name][int(query_align_len):])
                seen[ref_name]["start"].append([len(extra_start_seq), query_name, extra_start_seq, "end", "1", ref_name, percent_id, combined_query_name, query_to_combined_strand])
                report_list.append("Grabbing sequence in process_query_combined. 4th" + " " + query_name + " " + ref_name + " " + combined_query_name)
    elif strand == "-1" and (int(ref_start) < 24 or int(ref_end) > int(ref_len) - 23):
        end_dif = int(query_end)
        start_dif = int(query_len) - int(query_start)
        
        if end_dif > ref_end_dif:
            if query_overlap_end == "start" and query_to_combined_strand == "1":
                try:
                    extra_end_seq = fastaIO.reverse_complement(assembly[combined_query_name][:-(end_dif - ref_end_dif)])
                except:
                    print "Error with assembly. Info:\n", "\t".join([ref_start, ref_end, query_start, query_end, ref_align_len, query_align_len, percent_id, ref_len, query_len, ref_coverage, query_coverage, frame, strand, ref_name, query_name, tag])
                    raise
                seen[ref_name]["end"].append([len(extra_end_seq), query_name, extra_end_seq, "end", "-1", ref_name, percent_id, combined_query_name, query_to_combined_strand])
                report_list.append("Grabbing sequence in process_query_combined. 5th" + " " + query_name + " " + ref_name + " " + combined_query_name)
            elif query_overlap_end == "start" and query_to_combined_strand == "-1":
                extra_end_seq = assembly[combined_query_name][(end_dif - ref_end_dif):]
                seen[ref_name]["end"].append([len(extra_end_seq), query_name, extra_end_seq, "end", "-1", ref_name, percent_id, combined_query_name, query_to_combined_strand])
                report_list.append("Grabbing sequence in process_query_combined. 6th" + " " + query_name + " " + ref_name + " " + combined_query_name)
        if start_dif > ref_start_dif:
            if query_overlap_end == "end" and query_to_combined_strand == "1":
                extra_start_seq = fastaIO.reverse_complement(assembly[combined_query_name][int(query_align_len):])
                seen[ref_name]["start"].append([len(extra_start_seq), query_name, extra_start_seq, "start", "-1", ref_name, percent_id, combined_query_name, query_to_combined_strand])
                report_list.append("Grabbing sequence in process_query_combined. 7th" + " " + query_name + " " + ref_name + " " + combined_query_name)
            elif query_overlap_end == "end" and query_to_combined_strand == "-1":
                try:
                    extra_start_seq = assembly[combined_query_name][:-int(query_align_len)]
                except:
                    print "Combined query missing from assembly. ref_name and query_name :", query_name, ref_name
                    raise
                seen[ref_name]["start"].append([len(extra_start_seq), query_name, extra_start_seq, "start", "-1", ref_name, percent_id, combined_query_name, query_to_combined_strand])
                report_list.append("Grabbing sequence in process_query_combined. 8th" + " " + query_name + " " + ref_name + " " + combined_query_name)
    return bad, seen, assembly, good

def process_both_combined(seen, good, assembly, bad, report_list, contigs, ref_start, ref_end, query_start, query_end, ref_align_len, query_align_len, percent_id, ref_len, query_len, ref_coverage, query_coverage, frame, strand, ref_name, query_name, tag):
    if ref_name not in seen:
        seen[ref_name]["start"] = []
        seen[ref_name]["end"] = []
    
    combined_ref_name = good[ref_name][0]
    combined_ref_overlap_end = good[ref_name][1]
    ref_overlap_end = good[ref_name][2]
    ref_to_combined_strand = good[ref_name][3]
    
    combined_query_name = good[query_name][0]
    combined_query_overlap_end = good[query_name][1]
    query_overlap_end = good[query_name][2]
    query_to_combined_strand = good[query_name][3]
    
    ref_end_dif = int(ref_len) - int(ref_end)
    ref_start_dif = int(ref_start)
    if strand == "1" and (int(ref_start) < 24 or int(ref_end) > int(ref_len) - 23):
        end_dif = int(query_len) - int(query_end)
        start_dif = int(query_start)
        
        '''
        seen[ref_name][end of overlap relative to ref] = [len query extension, query_name, seq of extension, end of overlap relative to query, orientation of query to actual ref, actual ref_name to add extension to, which (ref or query) are combined, actual query_name, ori to actual query]
        '''
        
        if end_dif > ref_end_dif:
            if ref_overlap_end == "start" and query_overlap_end == "end" and ref_to_combined_strand == "1":
                if query_to_combined_strand == "1":
                    extra_end_seq = assembly[combined_query_name][(end_dif - ref_end_dif):]
                    seen[ref_name]["end"].append([len(extra_end_seq), query_name, extra_end_seq, "start", "1", combined_ref_name, percent_id, combined_query_name, query_to_combined_strand])
                    report_list.append("Grabbing sequence in process_both_combined. First" + " " + query_name + " " + ref_name + " " + combined_query_name + " " + combined_ref_name)                    
                else:
                    extra_end_seq = fastaIO.reverse_complement(assembly[combined_query_name][:-(end_dif - ref_end_dif)])
                    seen[ref_name]["end"].append([len(extra_end_seq), query_name, extra_end_seq, "start", "1", combined_ref_name, percent_id, combined_query_name, query_to_combined_strand])
                    report_list.append("Grabbing sequence in process_both_combined. Second" + " " + query_name + " " + ref_name + " " + combined_query_name + " " + combined_ref_name)
                    
            elif ref_overlap_end == "start" and query_overlap_end == "end" and ref_to_combined_strand == "-1":
                if query_to_combined_strand == "1":
                    extra_end_seq = fastaIO.reverse_complement(assembly[combined_query_name][:(end_dif - ref_end_dif)])
                    seen[ref_name]["start"].append([len(extra_end_seq), query_name, extra_end_seq, "end", "-1", combined_ref_name, percent_id, combined_query_name, query_to_combined_strand])
                    report_list.append("Grabbing sequence in process_both_combined. Third" + " " + query_name + " " + ref_name + " " + combined_query_name + " " + combined_ref_name)
                else:
                    extra_end_seq = assembly[combined_query_name][:-(end_dif - ref_end_dif)]
                    seen[ref_name]["start"].append([len(extra_end_seq), query_name, extra_end_seq, "end", "-1", combined_ref_name, percent_id, combined_query_name, query_to_combined_strand])
                    report_list.append("Grabbing sequence in process_both_combined. Fourth" + " " + query_name + " " + ref_name + " " + combined_query_name + " " + combined_ref_name)
        if start_dif > ref_start_dif:
            if ref_overlap_end == "end" and query_overlap_end == "start" and ref_to_combined_strand == "1":
                if query_to_combined_strand == "1":
                    extra_start_seq = assembly[combined_query_name][:-int(query_align_len)]
                    seen[ref_name]["start"].append([len(extra_start_seq), query_name, extra_start_seq, "start", "1", combined_ref_name, percent_id, combined_query_name, query_to_combined_strand])
                    report_list.append("Grabbing sequence in process_both_combined. Fifth" + " " + query_name + " " + ref_name + " " + combined_query_name + " " + combined_ref_name)
                else:
                    extra_start_seq = fastaIO.reverse_complement(assembly[combined_query_name][int(query_align_len):])
                    seen[ref_name]["start"].append([len(extra_start_seq), query_name, extra_start_seq, "start", "1", combined_ref_name, percent_id, combined_query_name, query_to_combined_strand])
                    report_list.append("Grabbing sequence in process_both_combined. Sixth" + " " + query_name + " " + ref_name + " " + combined_query_name + " " + combined_ref_name)
            elif ref_overlap_end == "end" and query_overlap_end == "start" and ref_to_combined_strand == "-1":
                if query_to_combined_strand == "1":
                    extra_start_seq = fastaIO.reverse_complement(assembly[combined_query_name][:-int(query_align_len)])
                    seen[ref_name]["end"].append([len(extra_start_seq), query_name, extra_start_seq, "start", "-1", combined_ref_name, percent_id, combined_query_name, query_to_combined_strand])
                    report_list.append("Grabbing sequence in process_both_combined. Seventh" + " " + query_name + " " + ref_name + " " + combined_query_name + " " + combined_ref_name)
                else:
                    extra_start_seq = assembly[combined_query_name][int(query_align_len):]
                    seen[ref_name]["end"].append([len(extra_start_seq), query_name, extra_start_seq, "start", "-1", combined_ref_name, percent_id, combined_query_name, query_to_combined_strand])
                    report_list.append("Grabbing sequence in process_both_combined. Eighth" + " " + query_name + " " + ref_name + " " + combined_query_name + " " + combined_ref_name)
    elif strand == "-1" and (int(ref_start) < 24 or int(ref_end) > int(ref_len) - 23):
        end_dif = int(query_end)
        start_dif = int(query_len) - int(query_start)
        
        '''
        seen[ref_name][end of overlap relative to ref] = [len query extension, query_name, seq of extension, end of overlap relative to query, orientation of query to actual ref, actual ref_name to add extension to, which (ref or query) are combined, actual query_name, ori to actual query]
        '''
        
        if end_dif > ref_end_dif:
            if ref_overlap_end == "start" and query_overlap_end == "start" and ref_to_combined_strand == "1":
                if query_to_combined_strand == "1":
                    extra_end_seq = fastaIO.reverse_complement(assembly[combined_query_name][:-(end_dif - ref_end_dif)])
                    seen[ref_name]["end"].append([len(extra_end_seq), query_name, extra_end_seq, "end", "-1", combined_ref_name, percent_id, combined_query_name, query_to_combined_strand])
                    report_list.append("Grabbing sequence in process_both_combined. 9th" + " " + query_name + " " + ref_name + " " + combined_query_name + " " + combined_ref_name)
                else:
                    extra_end_seq = assembly[combined_query_name][(end_dif - ref_end_dif):]
                    seen[ref_name]["end"].append([len(extra_end_seq), query_name, extra_end_seq, "end", "-1", combined_ref_name, percent_id, combined_query_name, query_to_combined_strand])
                    report_list.append("Grabbing sequence in process_both_combined. 10th" + " " + query_name + " " + ref_name + " " + combined_query_name + " " + combined_ref_name)
            elif ref_overlap_end == "start" and query_overlap_end == "start" and ref_to_combined_strand == "-1":
                if query_to_combined_strand == "1":
                    extra_end_seq = assembly[combined_query_name][:-(end_dif - ref_end_dif)]
                    seen[ref_name]["start"].append([len(extra_end_seq), query_name, extra_end_seq, "end", "1", combined_ref_name, percent_id, combined_query_name, query_to_combined_strand])
                    report_list.append("Grabbing sequence in process_both_combined. 11th" + " " + query_name + " " + ref_name + " " + combined_query_name + " " + combined_ref_name)
                else:
                    extra_end_seq = fastaIO.reverse_complement(assembly[combined_query_name][(end_dif - ref_end_dif):])
                    seen[ref_name]["start"].append([len(extra_end_seq), query_name, extra_end_seq, "end", "1", combined_ref_name, percent_id, combined_query_name, query_to_combined_strand])
                    report_list.append("Grabbing sequence in process_both_combined. 12th" + " " + query_name + " " + ref_name + " " + combined_query_name + " " + combined_ref_name)
                    
        if start_dif > ref_start_dif:
            if ref_overlap_end == "end" and query_overlap_end == "end" and ref_to_combined_strand == "1":
                if query_to_combined_strand == "1":
                    extra_start_seq = fastaIO.reverse_complement(assembly[combined_query_name][int(query_align_len):])
                    seen[ref_name]["start"].append([len(extra_start_seq), query_name, extra_start_seq, "start", "-1", combined_ref_name, percent_id, combined_query_name, query_to_combined_strand])
                    report_list.append("Grabbing sequence in process_both_combined. 13th" + " " + query_name + " " + ref_name + " " + combined_query_name + " " + combined_ref_name)
                else:
                    extra_start_seq = assembly[combined_query_name][:-int(query_align_len)]
                    seen[ref_name]["start"].append([len(extra_start_seq), query_name, extra_start_seq, "start", "-1", combined_ref_name, percent_id, combined_query_name, query_to_combined_strand])
                    report_list.append("Grabbing sequence in process_both_combined. 14th" + " " + query_name + " " + ref_name + " " + combined_query_name + " " + combined_ref_name)
            elif ref_overlap_end == "end" and query_overlap_end == "end" and ref_to_combined_strand == "-1":
                if query_to_combined_strand == "1":
                    extra_start_seq = fastaIO.reverse_complement(assembly[combined_query_name][int(query_align_len):])
                    seen[ref_name]["end"].append([len(extra_start_seq), query_name, extra_start_seq, "start", "1", combined_ref_name, percent_id, combined_query_name, query_to_combined_strand])
                    report_list.append("Grabbing sequence in process_both_combined. 15th" + " " + query_name + " " + ref_name + " " + combined_query_name + " " + combined_ref_name)
                else:
                    extra_start_seq = fastaIO.reverse_complement(assembly[combined_query_name][int(query_align_len):])
                    seen[ref_name]["end"].append([len(extra_start_seq), query_name, extra_start_seq, "start", "1", combined_ref_name, percent_id, combined_query_name, query_to_combined_strand])
                    report_list.append("Grabbing sequence in process_both_combined. 16th" + " " + query_name + " " + ref_name + " " + combined_query_name + " " + combined_ref_name)
    return bad, seen, assembly, good

def process_contigs(last_ref, final_name, final_name_ori, contig_self, contig_in, seen_entry):
    
    '''
    dict contig_in
    contig_in[base_contig_name][contained_contig_name] = ori of contained to base contig
    
    dict contig_self
    contig_self[contained_supercontgi_num][contained_contig_num][base_contig_name] = ori of contained to base contig
    '''
    
    for item in seen_entry:
        if item[5] == final_name:
            if "." in item[7]:
                super_num, contig_num = item[7].split(".")
                contig_self[super_num][contig_num][final_name] = final_name_ori
                contig_in[final_name][item[7]] = final_name_ori
                for item in contig_in[item[7]]:
                    if item not in contig_in[final_name]:
                        if (final_name_ori == "1" and contig_in[item[7]][item] == "1") or (final_name_ori == "-1" and contig_in[item[7]][item] == "-1"):
                            contig_in[final_name][item] = "1"
                        else:
                            contig_in[final_name][item] = "-1"
                contig_in.pop(item[7], None)
        elif item[7] == final_name:
            if "." in last_ref:
                super_num, contig_num = last_ref.split(".")
                contig_self[super_num][contig_num][final_name] = final_name_ori
                contig_in[final_name][last_ref] = contig_self[super_num][contig_num][final_name]
                for item in contig_in[item[5]]:
                    if item not in contig_in[final_name]:
                        if (final_name_ori == "1" and contig_in[item[5]][item] == "1") or (final_name_ori == "-1" and contig_in[item[5]][item] == "-1"):
                            contig_in[final_name][item] = "1"
                        else:
                            contig_in[final_name][item] = "-1"
                contig_in.pop(item[5], None)
        else:
            print "Not sure what is the final name comared to inputs.\nFinal_name =", final_name, "  last_ref =", last_ref, "  actual_ref =", item[5], "  query =", item[1], "  actual_query =", item[7]
            if "." in item[7]:
                super_num, contig_num = item[7].split(".")
                contig_self[super_num][contig_num][final_name] = final_name_ori
                contig_in[final_name][item[7]] = final_name_ori
                for item in contig_in[item[7]]:
                    if item not in contig_in[final_name]:
                        if (final_name_ori == "1" and contig_in[item[7]][item] == "1") or (final_name_ori == "-1" and contig_in[item[7]][item] == "-1"):
                            contig_in[final_name][item] = "1"
                        else:
                            contig_in[final_name][item] = "-1"
                contig_in.pop(item[7], None)
    
        for key in contig_in[final_name]:
            super_num, contig_num = key.split(".")
            screen_list = contig_self[super_num][contig_num].keys()
            for item in screen_list:
                if item in contig_in[final_name]:
                    contig_self[super_num][contig_num].pop(item, None)
    
    return contig_self, contig_in


def find_longest_extension(seen, good, bad, report_list, assembly, last_ref, contigs, covers, processed, contig_self, contig_in):
    final_name_list = []
    final_name_dict = {}
    if last_ref in good:
        actual_ref = good[last_ref][0]
    elif last_ref in covers:
        actual_ref = covers[last_ref]
        covers.pop(last_ref, None)
    else:
        actual_ref = last_ref
    if actual_ref in assembly and actual_ref not in bad:
        if actual_ref != last_ref:
            final_name_dict[actual_ref] = good[last_ref][3]
        else:
            final_name_dict[actual_ref] = "1"
        final_name_list.append(actual_ref)
        start_seq = ''
        end_seq = ''
        pop = 0
        if last_ref in bad:
            print "Last ref", last_ref, "is in bad. Top"
        if actual_ref in bad:
            print "Actual ref", actual_ref, "is in bad. Top"
        
        #sort overlap lists
        if len(seen[last_ref]["start"]) > 0:            
            seen[last_ref]["start"].sort(key=itemgetter(0), reverse=True)
            while len(seen[last_ref]["start"]) > 0 and seen[last_ref]["start"][0][1] in bad:
                seen[last_ref]["start"].pop(0)
            if len(seen[last_ref]["start"]) > 0:
                start_seq = seen[last_ref]["start"][0][2]
                report_list.append(last_ref + " (" + actual_ref + ")"+ "\t" + seen[last_ref]["start"][0][1] + " part of " + seen[last_ref]["start"][0][7] + "\t" + "start_extended_" + str(seen[last_ref]["start"][0][0]))
                pop += 1 
                if seen[last_ref]["start"][0][5] not in final_name_dict:
                    final_name_dict[seen[last_ref]["start"][0][5]] = "1"
                    final_name_list.append(seen[last_ref]["start"][0][5])
                if seen[last_ref]["start"][0][7] not in final_name_dict:    
                    if seen[last_ref]["start"][0][4] == "1":
                        if seen[last_ref]["start"][0][8] == "1":
                            final_name_dict[seen[last_ref]["start"][0][7]] = "1"
                        else:
                            final_name_dict[seen[last_ref]["start"][0][7]] = "-1"
                    else:
                        if seen[last_ref]["start"][0][8] == "-1":
                            final_name_dict[seen[last_ref]["start"][0][7]] = "1"
                        else:
                            final_name_dict[seen[last_ref]["start"][0][7]] = "-1"
                    final_name_list.append(seen[last_ref]["start"][0][7])
        
        if len(seen[last_ref]["end"]) > 0:            
            seen[last_ref]["end"].sort(key=itemgetter(0), reverse=True)
            while len(seen[last_ref]["end"]) > 0 and seen[last_ref]["end"][0][1] in bad:
                seen[last_ref]["end"].pop(0)
            if len(seen[last_ref]["end"]) > 0:
                end_seq = seen[last_ref]["end"][0][2]
                report_list.append(last_ref + " (" + actual_ref + ")"+ "\t" + seen[last_ref]["end"][0][1] + " part of " + seen[last_ref]["end"][0][7] + "\t" + "end_extended_" + str(seen[last_ref]["end"][0][0]))
                pop += 2
                if seen[last_ref]["end"][0][7] not in final_name_dict:
                    if seen[last_ref]["end"][0][4] == "1":
                        if seen[last_ref]["end"][0][8] == "1":
                            final_name_dict[seen[last_ref]["end"][0][7]] = "1"
                        else:
                            final_name_dict[seen[last_ref]["end"][0][7]] = "-1"
                    else:
                        if seen[last_ref]["end"][0][8] == "-1":
                            final_name_dict[seen[last_ref]["end"][0][7]] = "1"
                        else:
                            final_name_dict[seen[last_ref]["end"][0][7]] = "-1"
                    final_name_list.append(seen[last_ref]["end"][0][7])
            
        if pop == 0:
            seen = Vividict()
            last_ref = ''
            return seen, assembly, good, bad, last_ref, contigs, covers, processed
        final_name_list.sort()
        final_name = final_name_list.pop(0)
        
        '''        
    seen dictionary:
    seen[ref_name][end of overlap relative to ref] = [len query extension, query_name, seq of extension, end of overlap relative to query, orientation of query to actual ref, actual ref_name to add extension to, which (ref or query) are combined, actual query_name, ori to actual query]
    
    good dictionary:
    good[best_query_contig] = [ref_contig, end of ref_contig extended, end of query_contig used in extension, orientation of query to final orientation (either just to ref if it's a single contig or to the orientation of merged contigs)]
        '''
            
        if pop == 3:
            if seen[last_ref]["start"][0][1] == seen[last_ref]["end"][0][1] and seen[last_ref]["start"][0][7] == seen[last_ref]["start"][0][1] and seen[last_ref]["end"][0][7] == seen[last_ref]["end"][0][1]:
                if last_ref != seen[last_ref]["start"][0][5]:
                    print "Problem: last_ref is not equal to actual_ref in seen even though the query covers the ref. Last_ref=", last_ref, "actual_ref in seen=", seen[last_ref]["start"][0][5] 
                covers[seen[last_ref]["start"][0][1]] = actual_ref
                report_list.append(last_ref + " (" + actual_ref + ")"+ "\t" + seen[last_ref]["start"][0][1] + "\t" + "query covers ref")
                if "." in seen[last_ref]["start"][0][1]:
                    super_num, contig_num = seen[last_ref]["start"][0][1].split(".")
                    contig_self[super_num][contig_num] = {[last_ref]: seen[last_ref]["start"][4]}
                    contig_in[last_ref][seen[last_ref]["start"][0][1]] = seen[last_ref]["start"][4]
                contigs[actual_ref] = contigs[seen[last_ref]["start"][0][1]]
                contigs.pop(seen[last_ref]["start"][0][1], None)
                contigs.pop(seen[last_ref]["start"][0][7], None)
                assembly.pop(seen[last_ref]["start"][0][1], None)
                contigs.pop(seen[last_ref]["end"][0][1], None)
                contigs.pop(seen[last_ref]["end"][0][7], None)
                processed[last_ref] = 1

            elif seen[last_ref]["start"][0][7] == seen[last_ref]["end"][0][7]:
                report_list.append(last_ref + " (" + actual_ref + ")"+ "\t" + seen[last_ref]["start"][0][1] + " and " + seen[last_ref]["start"][0][1] + " both in " + seen[last_ref]["start"][0][7] + "\t" + " Combined query covers ref")
                start_seq = ''
                end_seq = ''
                assembly[actual_ref] = assembly[seen[last_ref]["start"][0][7]]
                contigs[actual_ref] = contigs[seen[last_ref]["start"][0][7]]
                if actual_ref < seen[last_ref]["start"][0][7]:
                    covers[seen[last_ref]["start"][0][7]] = actual_ref
                    assembly.pop(seen[last_ref]["start"][0][7], None)
                    assembly.pop(seen[last_ref]["start"][0][1], None)
                    assembly.pop(seen[last_ref]["end"][0][1], None)
                    contigs.pop(seen[last_ref]["start"][0][7], None)
                    contigs.pop(seen[last_ref]["start"][0][1], None)
                    contigs.pop(seen[last_ref]["end"][0][1], None)
                contig_self, contig_in = process_contigs(last_ref, final_name, final_name_dict[final_name], contig_self, contig_in, seen[last_ref]["start"])
                processed[last_ref] = 1
                
            else:
                contig_self, contig_in = process_contigs(last_ref, final_name, final_name_dict[final_name], contig_self, contig_in, seen[last_ref]["start"])
                if seen[last_ref]["start"][0][1] in good:
                    if seen[last_ref]["start"][0][1] != actual_ref:
                        bad[seen[last_ref]["start"][0][1]] = 1
                        good.pop(seen[last_ref]["start"][0][1], None)
                        assembly.pop(seen[last_ref]["start"][0][1], None)
                        if seen[last_ref]["start"][0][1] != seen[last_ref]["start"][0][7]:
                            if (seen[last_ref]["start"][0][4] == "1" and seen[last_ref]["start"][0][8] == "-1") or (seen[last_ref]["start"][0][4] == "-1" and seen[last_ref]["start"][0][8] == "1"):
                                contigs[seen[last_ref]["start"][0][7]].reverse()
                        elif seen[last_ref]["start"][0][4] == "-1":
                            contigs[seen[last_ref]["start"][0][7]].reverse()
                        contigs[seen[last_ref]["start"][0][7]].extend(contigs[actual_ref])
                        contigs[actual_ref] = contigs[seen[last_ref]["start"][0][7]]
                        if last_ref in good:
                            good.pop(last_ref, None)
                        contigs.pop(seen[last_ref]["start"][0][1], None)
                        if seen[last_ref]["start"][0][7] != actual_ref:
                            assembly.pop(seen[last_ref]["start"][0][7], None)
                            contigs.pop(seen[last_ref]["start"][0][7], None)

                else:
                    good[seen[last_ref]["start"][0][1]] = [actual_ref, "start", seen[last_ref]["start"][0][3], seen[last_ref]["start"][0][4]]
                    if seen[last_ref]["start"][0][1] != seen[last_ref]["start"][0][7]:
                        if (seen[last_ref]["start"][0][4] == "1" and seen[last_ref]["start"][0][8] == "-1") or (seen[last_ref]["start"][0][4] == "-1" and seen[last_ref]["start"][0][8] == "1"):
                            contigs[seen[last_ref]["start"][0][7]].reverse()
                    elif seen[last_ref]["start"][0][4] == "-1":
                        contigs[seen[last_ref]["start"][0][7]].reverse()
                    try:
                        contigs[seen[last_ref]["start"][0][7]].extend(contigs[actual_ref])
                    except:
                        print "last_ref: ", last_ref, "actual_ref:", actual_ref, "\nseen[last_ref]['start'] =", seen[last_ref]["start"], "\ngood[seen[last_ref]['start'][0][1]] =", good[seen[last_ref]["start"][0][1]], "\nbad[seen[last_ref]['start'][0][7]] = ", bad[seen[last_ref]["start"][0][7]]
                        raise
                    contigs[actual_ref] = contigs[seen[last_ref]["start"][0][7]]
                    if seen[last_ref]["start"][0][1] != actual_ref:
                        contigs.pop(seen[last_ref]["start"][0][1], None)
                        assembly.pop(seen[last_ref]["start"][0][1], None)
                    if seen[last_ref]["start"][0][7] != actual_ref:
                        assembly.pop(seen[last_ref]["start"][0][7], None)
                        contigs.pop(seen[last_ref]["start"][0][7], None)
                    if last_ref in good:
                        good.pop(last_ref, None)
                
                contig_self, contig_in = process_contigs(last_ref, final_name, final_name_dict[final_name], contig_self, contig_in, seen[last_ref]["end"])
                if seen[last_ref]["end"][0][1] in good:
                    if seen[last_ref]["end"][0][1] != actual_ref:
                        bad[seen[last_ref]["end"][0][1]] = 1
                        good.pop(seen[last_ref]["end"][0][1], None)
                        assembly.pop(seen[last_ref]["end"][0][1], None)
                        if seen[last_ref]["end"][0][1] != seen[last_ref]["end"][0][7]:
                            if (seen[last_ref]["end"][0][4] == "1" and seen[last_ref]["end"][0][8] == "-1") or (seen[last_ref]["end"][0][4] == "-1" and seen[last_ref]["end"][0][8] == "1"):
                                contigs[seen[last_ref]["end"][0][7]].reverse()
                        elif seen[last_ref]["end"][0][4] == "-1":
                            contigs[seen[last_ref]["end"][0][7]].reverse()
                        try:
                            contigs[actual_ref].extend(contigs[seen[last_ref]["end"][0][7]])
                        except:
                            print "last_ref: ", last_ref, "actual_ref:", actual_ref, "\nseen[last_ref] =", seen[last_ref]["end"]
                            raise    
                            
                        if last_ref in good:
                            good.pop(last_ref, None)
                        contigs.pop(seen[last_ref]["end"][0][1], None)
                        if seen[last_ref]["end"][0][7] != actual_ref:
                            contigs.pop(seen[last_ref]["end"][0][7], None)
                            assembly.pop(seen[last_ref]["end"][0][7], None)

                else:
                    good[seen[last_ref]["end"][0][1]] = [actual_ref, "end", seen[last_ref]["end"][0][3], seen[last_ref]["end"][0][4]]
                    #report_list.append(last_ref + " (" + actual_ref + ")"+ "\t" + seen[last_ref]["end"][0][1] + "\t" + "query added to ref end")
                    if seen[last_ref]["end"][0][1] != seen[last_ref]["end"][0][7]:
                        if (seen[last_ref]["end"][0][4] == "1" and seen[last_ref]["end"][0][8] == "-1") or (seen[last_ref]["end"][0][4] == "-1" and seen[last_ref]["end"][0][8] == "1"):
                            contigs[seen[last_ref]["end"][0][7]].reverse()
                    elif seen[last_ref]["end"][0][4] == "-1":
                        contigs[seen[last_ref]["end"][0][7]].reverse()
                    contigs[actual_ref].extend(contigs[seen[last_ref]["end"][0][7]])
                    if seen[last_ref]["end"][0][1] != actual_ref:
                        contigs.pop(seen[last_ref]["end"][0][1], None)
                        assembly.pop(seen[last_ref]["end"][0][1], None)
                    if seen[last_ref]["end"][0][7] != actual_ref:
                        assembly.pop(seen[last_ref]["end"][0][7], None)
                        contigs.pop(seen[last_ref]["end"][0][7], None)
                    if last_ref in good:
                        good.pop(last_ref, None)
                
            '''        
        seen dictionary:
        seen[ref_name][end of overlap relative to ref] = [len query extension, query_name, seq of extension, end of overlap relative to query, orientation of query to actual ref, actual ref_name to add extension to, which (ref or query) are combined, actual query_name, ori to actual query]
        
        good dictionary:
        good[best_query_contig] = [ref_contig, end of ref_contig extended, end of query_contig used in extension, orientation of query to final orientation (either just to ref if it's a single contig or to the orientation of merged contigs)]
            '''
        
        elif pop == 1:
            contig_self, contig_in = process_contigs(last_ref, final_name, final_name_dict[final_name], contig_self, contig_in, seen[last_ref]["start"])
            if seen[last_ref]["start"][0][1] in good:
                if seen[last_ref]["start"][0][1] != actual_ref:
                    bad[seen[last_ref]["start"][0][1]] = 1
                    good.pop(seen[last_ref]["start"][0][1], None)
                    assembly.pop(seen[last_ref]["start"][0][1], None)
                    if seen[last_ref]["start"][0][1] != seen[last_ref]["start"][0][7]:
                        if (seen[last_ref]["start"][0][4] == "1" and seen[last_ref]["start"][0][8] == "-1") or (seen[last_ref]["start"][0][4] == "-1" and seen[last_ref]["start"][0][8] == "1"):
                            contigs[seen[last_ref]["start"][0][7]].reverse()
                    elif seen[last_ref]["start"][0][4] == "-1":
                        contigs[seen[last_ref]["start"][0][7]].reverse()
                    contigs[seen[last_ref]["start"][0][7]].extend(contigs[actual_ref])
                    contigs[actual_ref] = contigs[seen[last_ref]["start"][0][7]]
                    if last_ref in good:
                        good.pop(last_ref, None)
                    contigs.pop(seen[last_ref]["start"][0][1], None)
                    if seen[last_ref]["start"][0][7] != actual_ref:
                        assembly.pop(seen[last_ref]["start"][0][7], None)
                        contigs.pop(seen[last_ref]["start"][0][7], None)
            else:
                good[seen[last_ref]["start"][0][1]] = [actual_ref, "start", seen[last_ref]["start"][0][3], seen[last_ref]["start"][0][4]]
                if seen[last_ref]["start"][0][1] != seen[last_ref]["start"][0][7]:
                    if (seen[last_ref]["start"][0][4] == "1" and seen[last_ref]["start"][0][8] == "-1") or (seen[last_ref]["start"][0][4] == "-1" and seen[last_ref]["start"][0][8] == "1"):
                        contigs[seen[last_ref]["start"][0][7]].reverse()
                elif seen[last_ref]["start"][0][4] == "-1":
                    contigs[seen[last_ref]["start"][0][7]].reverse()
                contigs[seen[last_ref]["start"][0][7]].extend(contigs[actual_ref])
                contigs[actual_ref] = contigs[seen[last_ref]["start"][0][7]]
                if seen[last_ref]["start"][0][1] != actual_ref:
                    contigs.pop(seen[last_ref]["start"][0][1], None)
                    assembly.pop(seen[last_ref]["start"][0][1], None)
                if seen[last_ref]["start"][0][7] != actual_ref:
                    assembly.pop(seen[last_ref]["start"][0][7], None)
                    contigs.pop(seen[last_ref]["start"][0][7], None)
                if last_ref in good:
                    good.pop(last_ref, None)
                        
        elif pop == 2:
            contig_self, contig_in = process_contigs(last_ref, final_name, final_name_dict[final_name], contig_self, contig_in, seen[last_ref]["end"])
            if seen[last_ref]["end"][0][1] in good:
                if seen[last_ref]["end"][0][1] != actual_ref:
                    bad[seen[last_ref]["end"][0][1]] = 1
                    if seen[last_ref]["end"][0][1] != seen[last_ref]["end"][0][7]:
                        if (seen[last_ref]["end"][0][4] == "1" and seen[last_ref]["end"][0][8] == "-1") or (seen[last_ref]["end"][0][4] == "-1" and seen[last_ref]["end"][0][8] == "1"):
                            contigs[seen[last_ref]["end"][0][7]].reverse()
                    elif seen[last_ref]["end"][0][4] == "-1":
                        contigs[seen[last_ref]["end"][0][7]].reverse()
                    contigs[actual_ref].extend(contigs[seen[last_ref]["end"][0][7]])
                    good.pop(seen[last_ref]["end"][0][1], None)
                    assembly.pop(seen[last_ref]["end"][0][1], None)
                    contigs.pop(seen[last_ref]["end"][0][1], None)
                    if seen[last_ref]["end"][0][7] != actual_ref:
                        assembly.pop(seen[last_ref]["end"][0][7], None)
                        contigs.pop(seen[last_ref]["end"][0][7], None)
                    if last_ref in good:
                        good.pop(last_ref, None)

            else:
                good[seen[last_ref]["end"][0][1]] = [actual_ref, "end", seen[last_ref]["end"][0][3], seen[last_ref]["end"][0][4]]
                if seen[last_ref]["end"][0][1] != seen[last_ref]["end"][0][7]:
                    if (seen[last_ref]["end"][0][4] == "1" and seen[last_ref]["end"][0][8] == "-1") or (seen[last_ref]["end"][0][4] == "-1" and seen[last_ref]["end"][0][8] == "1"):
                        contigs[seen[last_ref]["end"][0][7]].reverse()
                elif seen[last_ref]["end"][0][4] == "-1":
                    contigs[seen[last_ref]["end"][0][7]].reverse()
                contigs[actual_ref].extend(contigs[seen[last_ref]["end"][0][7]])
                if seen[last_ref]["end"][0][1] != actual_ref:
                    contigs.pop(seen[last_ref]["end"][0][1], None)
                    assembly.pop(seen[last_ref]["end"][0][1], None)
                if seen[last_ref]["end"][0][7] != actual_ref:
                    assembly.pop(seen[last_ref]["end"][0][7], None)
                    contigs.pop(seen[last_ref]["end"][0][7], None)
                if last_ref in good:
                    good.pop(last_ref, None)
        
        seen[last_ref]["start"].pop(0, None)
        seen[last_ref]["end"].pop(0, None)
        try:
            new_seq = start_seq + assembly[actual_ref] + end_seq
        except:
            print "actual_ref not in assebly dict! actual_ref =", actual_ref, "last_ref =", last_ref, "/nbad[actual_ref] =", bad[actual_ref], "\ngood[actual_ref] =", good[actual_ref]
            raise
        
        if final_name != actual_ref:
            bad[actual_ref] = 1
            contigs[final_name] = contigs[actual_ref]
            report_list.append(actual_ref + "\t" + final_name + "\t" + "actual_ref is not equal to final_ref")
            if contigs[final_name][0] in good:
                good[contigs[final_name][0]][0] = final_name
            if contigs[final_name][-1] in good:
                good[contigs[final_name][-1]][0] = final_name
            contigs.pop(actual_ref, None)
            assembly.pop(actual_ref, None)
        if final_name_dict[final_name] == "1":
            assembly[final_name] = new_seq
        else:
            assembly[final_name] = fastaIO.reverse_complement(new_seq)
            contigs[final_name].reverse()
        processed[last_ref] = 1    
            
        for item in seen[last_ref]["start"]:
            if item[1] not in bad:
                if item[1] not in good:
                    both = 0
                    bad[item[1]] = 1
                    assembly.pop(item[1], None)
                    contigs.pop(item[1], None)
                    for item2 in seen[last_ref]["end"]:
                        if item2[1] == item[1]:
                            report_list.append(last_ref + "\t" + item[1] + "\t" + "both_ends_short")
                            both = 1
                    if both == 0:
                        report_list.append(last_ref + "\t" + item[1] + "\t" + "start_short")
                else:
                    report_list.append(last_ref + "\t" + item[1] + "\t" + "start_short but contig is in good")
        
        for item in seen[last_ref]["end"]:
            if item[1] not in bad:
                if item[1] not in good:
                    bad[item[1]] = 1
                    assembly.pop(item[1], None)
                    contigs.pop(item[1], None)
                    report_list.append(last_ref + "\t" + item[1] + "\t" + "end_short")
                else:
                    report_list.append(last_ref + "\t" + item[1] + "\t" + "end_short but contig is in good")
                    
        for seq_name in final_name_list:
            if seq_name != final_name and seq_name not in good:
                bad[seq_name] = 1
                assembly.pop(seq_name, None)
                good.pop(seq_name, None)
                contigs.pop(seq_name, None)
        
    seen = Vividict()
    last_ref = ''
    return seen, assembly, good, bad, last_ref, contigs, covers, processed, contig_self, contig_in

def clear_multiple_matches(previous_ref, seen, bad, processed):
    last_ref = previous_ref
    if len(seen[last_ref]["start"]) > 0:
        if len(seen[last_ref]["end"]) > 0:
            if seen[last_ref]["start"][-1][1] == seen[last_ref]["end"][-1][1]:
                bad[last_ref] = 1
                processed[last_ref] = 1
                bad[seen[last_ref]["start"][-1][1]] = 1
                seen[last_ref]["start"].pop(-1)
                seen[last_ref]["end"].pop(-1)
            elif len(seen[last_ref]["start"]) > 1 and seen[last_ref]["start"][-1][1] == seen[last_ref]["start"][-2][1]:
                if seen[last_ref]["start"][-1][6] > seen[last_ref]["start"][-2][6]:
                    seen[last_ref]["start"].pop(-2)
                elif seen[last_ref]["start"][-1][6] == seen[last_ref]["start"][-2][6]:
                    if seen[last_ref]["start"][-1][0] >= seen[last_ref]["start"][-2][0]:
                        seen[last_ref]["start"].pop(-2)
                    else:
                        seen[last_ref]["start"].pop(-1)
                else:
                    seen[last_ref]["start"].pop(-1)
            elif len(seen[last_ref]["end"]) > 1 and seen[last_ref]["end"][-1][1] == seen[last_ref]["end"][-2][1]:
                if seen[last_ref]["end"][-1][6] > seen[last_ref]["end"][-2][6]:
                    seen[last_ref]["end"].pop(-2)
                elif seen[last_ref]["end"][-1][6] == seen[last_ref]["end"][-2][6]:
                    if seen[last_ref]["end"][-1][0] >= seen[last_ref]["end"][-2][0]:
                        seen[last_ref]["end"].pop(-2)
                    else:
                        seen[last_ref]["end"].pop(-1)
                else:
                    seen[last_ref]["end"].pop(-1)
        elif len(seen[last_ref]["start"]) > 1 and seen[last_ref]["start"][-1][1] == seen[last_ref]["start"][-2][1]:
            if seen[last_ref]["start"][-1][6] > seen[last_ref]["start"][-2][6]:
                seen[last_ref]["start"].pop(-2)
            elif seen[last_ref]["start"][-1][6] == seen[last_ref]["start"][-2][6]:
                if seen[last_ref]["start"][-1][0] >= seen[last_ref]["start"][-2][0]:
                    seen[last_ref]["start"].pop(-2)
                else:
                    seen[last_ref]["start"].pop(-1)
            else:
                seen[last_ref]["start"].pop(-1)
    elif len(seen[last_ref]["end"]) > 1:
        if seen[last_ref]["end"][-1][1] == seen[last_ref]["end"][-2][1]:
            if seen[last_ref]["end"][-1][6] > seen[last_ref]["end"][-2][6]:
                seen[last_ref]["end"].pop(-2)
            elif seen[last_ref]["end"][-1][6] == seen[last_ref]["end"][-2][6]:
                if seen[last_ref]["end"][-1][0] >= seen[last_ref]["end"][-2][0]:
                    seen[last_ref]["end"].pop(-2)
                else:
                    seen[last_ref]["end"].pop(-1)
            else:
                seen[last_ref]["end"].pop(-1)
    return seen, bad, processed

def supercontiger(contigs, assembly, contig_in, contig_self):
    
    '''
    dict contig_in
    contig_in[base_contig_name][contained_contig_name] = ori of contained to base contig
    
    dict contig_self
    contig_self[contained_supercontgi_num][contained_contig_num][base_contig_name] = ori of contained to base contig
    
    Ordered_dict contigs
    contigs[base_contig] = [ordered list of contigs contributing sequence to base contig]
    
    '''
    
    connection_track = {}
    
    for contig in contigs:
        if len(contigs[contig]) > 1:
            pass
    
def main():

    bad = {}
    good = {}
    processed = {}
    covers = {}
    report_list = []
    assembly = OrderedDict()
    contigs = OrderedDict()
    contig_self = Vividict()
    contig_track = []
    contig_in = Vividict()
    last = []
    seen = Vividict()
    bad_out = os.path.splitext(sys.argv[1])[0] + "_" + sys.argv[3] + "_report.out"
    contigs_out = os.path.splitext(sys.argv[1])[0] + "_" + sys.argv[3] + "_contigs.out"
    contig_in_out = os.path.splitext(sys.argv[1])[0] + "_" + sys.argv[3] + "_contig_in.out"
    contig_self_out = os.path.splitext(sys.argv[1])[0] + "_" + sys.argv[3] + "_contig_self.out"
    
    with open(sys.argv[2], "r") as f:
        for title, seq in fastaIO.FastaGeneralIterator(f):
            assembly[title] = seq
            contigs[title] = [title]
    
    with open(sys.argv[1], "r") as f:
        c = 0
        last_ref = ''
        last_query = ''
        for line in f:
            line = line.strip()
            if c == 0:
                c += 1
                continue
            if "CONTAINS" not in line and "IDENTITY" not in line and "END" not in line and "BEGIN" not in line and "CONTAINED" not in line:
                continue
            
            ref_start, ref_end, query_start, query_end, ref_align_len, query_align_len, percent_id, ref_len, query_len, ref_coverage, query_coverage, frame, strand, ref_name, query_name, tag = line.split("\t")
                        
            if last_ref and ref_name != last_ref:
                seen, assembly, good, bad, last_ref, contigs, covers, processed, contig_self, contig_in = find_longest_extension(seen, good, bad, report_list, assembly, last_ref, contigs, covers, processed, contig_self, contig_in)
            
            if ref_name == query_name:
                continue
            if ref_name in bad:
                continue
            if query_name in bad or query_name in covers or query_name in processed:
                continue
            if float(percent_id) <= 94.99:
                continue
            previous_ref = last_ref
            last_ref = ref_name
                        
            if ref_name not in good:
                if query_name not in good:
                    if query_name < ref_name:
                        continue
                    bad, seen, assembly, good, contigs, contig_self, contig_track = process_single(seen, good, assembly, bad, report_list, contigs, contig_self, contig_track, ref_start, ref_end, query_start, query_end, ref_align_len, query_align_len, percent_id, ref_len, query_len, ref_coverage, query_coverage, frame, strand, ref_name, query_name, tag)
                    if last_query == query_name and previous_ref == ref_name:
                        if len(seen[ref_name]["start"]) > 0 or len(seen[ref_name]["end"])> 0:
                            seen, bad, processed = clear_multiple_matches(previous_ref, seen, bad, processed)
                    last_query = query_name
                else:
                    if good[query_name][0] in bad:
                        continue
                    bad, seen, assembly, good = process_query_combined(seen, good, assembly, bad, report_list, contigs, ref_start, ref_end, query_start, query_end, ref_align_len, query_align_len, percent_id, ref_len, query_len, ref_coverage, query_coverage, frame, strand, ref_name, query_name, tag)
                    if last_query == query_name and previous_ref == ref_name:
                        if len(seen[ref_name]["start"]) > 0 or len(seen[ref_name]["end"])> 0:
                            seen, bad, processed = clear_multiple_matches(previous_ref, seen, bad, processed)
                    last_query = query_name
                            
            else: #ref_name in good:
                if good[ref_name][0] in bad or good[ref_name][0] == query_name:
                    continue
                if query_name in good:
                    if good[query_name][0] in bad:
                         continue
                    bad, seen, assembly, good = process_both_combined(seen, good, assembly, bad, report_list, contigs, ref_start, ref_end, query_start, query_end, ref_align_len, query_align_len, percent_id, ref_len, query_len, ref_coverage, query_coverage, frame, strand, ref_name, query_name, tag)
                    if last_query == query_name and previous_ref == ref_name:
                        if len(seen[ref_name]["start"]) > 0 or len(seen[ref_name]["end"])> 0:
                            seen, bad, processed = clear_multiple_matches(previous_ref, seen, bad, processed)
                    last_query = query_name
                else:
                    if query_name in bad:
                        continue
                    else:
                        bad, seen, assembly, good, contig_self, contig_track = process_ref_combined(seen, good, assembly, bad, report_list, contigs, contig_self, contig_track, ref_start, ref_end, query_start, query_end, ref_align_len, query_align_len, percent_id, ref_len, query_len, ref_coverage, query_coverage, frame, strand, ref_name, query_name, tag)
                        if last_query == query_name and previous_ref == ref_name:
                            if len(seen[ref_name]["start"]) > 0 or len(seen[ref_name]["end"])> 0:
                                seen, bad, processed = clear_multiple_matches(previous_ref, seen, bad, processed)
                        last_query = query_name
        
        seen, assembly, good, bad, last_ref, contigs, covers, processed, contig_self, contig_in = find_longest_extension(seen, good, bad, report_list, assembly, last_ref, contigs, covers, processed, contig_self, contig_in)
    if sys.argv[4].upper()[0] == "Y":
        contigs, assembly = supercontiger(contigs, assembly, contig_in, contig_self)
    else:
        with open(contig_in_out, "w", 1) as out:
            for contig in contig_in:
                print>>out, contig + "\t" + ",".join(contig_in[contig])
        with open(contig_self_out, "w", 1) as out:
            for super_num in contig_self:
                for contig_num in contig_self[super_num]:
                    print>>out, super_num + "." + contig_num + "\t" + ",".join(contig_self[super_num][contig_num])
        
    for item in good:
        if item not in contigs:
            report_list.append(item + "\t" + "Still in good but not contigs after processing\n")
        else:
            report_list.append(item + "\t" + "Still in good and contigs after processing\n")
        report_list.append(good[item])
        assembly.pop(item, None)
        contigs.pop(item, None)
        
    with open(bad_out, "w", 1) as out:
        for item in report_list:
            print>>out, item
    report_list = []
    
    assembly_out = os.path.splitext(sys.argv[2])[0] + "_" + sys.argv[3] + ".fa"
    with open(assembly_out, "w", 1) as out:
        for title in assembly:
            print>>out, ">" + title + "\n" + assembly[title]
    assembly = {}
    
    newdct = []
    for key in contigs:
        newdct.append(key)
    newdct.sort()       
    with open(contigs_out, "w", 1) as out:
        for item in newdct:
            print>>out, item + "\t" + "\t".join(contigs[item])
    
	return 0

if __name__ == '__main__':
	main()
