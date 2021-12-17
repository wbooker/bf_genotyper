import pysam
import os
import pathlib
from itertools import chain
import random

def calc_query_pos_from_cigar(cigar, strand):
    """
    copied from samplot source code (https://github.com/ryanlayer/samplot)
    Uses the CIGAR string to determine the query position of a read

    The cigar arg is a string like the following: 86M65S
    The strand arg is a boolean, True for forward strand and False for
    reverse

    Returns pair of ints for query start, end positions
    """

    cigar_ops = [[int(op[0]), op[1]] for op in re.findall("(\d+)([A-Za-z])", cigar)]

    order_ops = cigar_ops
    if not strand:  # - strand
        order_ops = order_ops[::-1]

    qs_pos = 0
    qe_pos = 0
    q_len = 0

    for op_position in range(len(cigar_ops)):
        op_len = cigar_ops[op_position][0]
        op_type = cigar_ops[op_position][1]

        if op_position == 0 and (op_type == "H" or op_type == "S"):
            qs_pos += op_len
            qe_pos += op_len
            q_len += op_len
        elif op_type == "H" or op_type == "S":
            q_len += op_len
        elif op_type == "M" or op_type == "I" or op_type == "X":
            qe_pos += op_len
            q_len += op_len

    return qs_pos, qe_pos


def get_split_event_type(primary, supplemental):
    """
        modified from the samplot source code (https://github.com/ryanlayer/samplot)

        Decide what type of event the read supports (del/normal, dup, inv)
        
        primary = [primary_reference, primary_strand, primary_q_start, primary_start_pos, primary_end_pos, bam_file.get_tid(primary_reference)]
        supplemental = [supp_reference, supp_strand, supp_q_start, supp_start_pos, supp_end_pos, bam_file.get_tid(supp_reference)]

        Here we first check positions and whether the read that starts second is completely overlapped by the left-most starting read,
        then we check the strand of each read. Then, we check that the query postion start, which tells us which read has trailing soft/hard 
        clipped bases and is therefore the read mapped after the initial read hit the breakpoint for forward reads, opposite for reverse.
        Then we check if first end extends past second read start, and if that with the query position matches duplication or deletion.
    """
    ### in the original samplot, a table with the primary and supplemental reads are ordered by start position before
    ### being input into this function, so we do this at the start to accurately check if one read is completely overlapped
    ### by another read. 
    if primary[3] < supplemental[3]:
        first = primary
        second = supplemental
    else:
        first = supplemental
        second = primary
    ### in samplot this checks first_pre_start against second_pre_end, but I don't see how this tells us anything and led
    ### to bad calls. modified checks to look at 1st read end overlapping 2nd read start

    event_type_by_strand_and_order = {
        (True, False): "INV",  # mixed strands
        (False, True): "INV",  # mixed strands
        (True, True, True, False): "DEL",  # forward strand
        (True, True, True, True): "DUP",  # forward strand
        (True, True, False, False): "DUP",  # forward strand
        (True, True, False, True): "DUP",  # forward strand
        (False, False, False, False): "DUP",  # reverse strand
        (False, False, False, True): "DUP",  # reverse strand
        (False, False, True, True): "DUP",  # reverse strand
        (False, False, True, False): "DEL",  # reverse strand
    }
    ## append first and second strand to orientations
    orientations = [first[1], second[1]]

    # if same strand, need query position info and whether second read overlaps first completely
    if orientations[0] == orientations[1]:

        # first query position smaller than second query position
        orientations.append(first[2] < second[2])

        # need to check if first read end extends past second read start (overlapping)
        orientations.append(first[4] > second[3])
    event_type = event_type_by_strand_and_order[tuple(orientations)]

    if first[0] != second[0]:
        if event_type == "INV":
            event_type = "InterChrmInversion"
        else:
            event_type = "InterChrm"

    return event_type


def get_breakpoints(primary, supplemental, split_type):
    """
    identify SV breakpoints using position and query_start position
    
    primary = [primary_reference, primary_strand, primary_q_start, primary_start_pos, primary_end_pos, bam_file.get_tid(primary_reference)]
    supplemental = [supp_reference, supp_strand, supp_q_start, supp_start_pos, supp_end_pos, bam_file.get_tid(supp_reference)]

    split_type = DUP, DEL, INV
    """
    if primary[3] < supplemental[3]:
        first = primary
        second = supplemental
    else:
        first = supplemental
        second = primary

    if split_type == "DUP":
        if first[2] < second[2]:
            left_bp = second[3]
            right_bp = first[4]
        else:
            left_bp = first[3]
            right_bp = second[4]
    return left_bp, right_bp 


def robust_get_MC(read):
    """
    sometimes MC cigar string missing for low mapping quality reads (mapq = 0). not really important as these are ignored, but 
    need to return something for functionality
    """  
    try:  
        return read.get_tag("MC")
    except KeyError:
        return str(random.randint(1,150))+"M"


def main():
    vcf_in = pysam.VariantFile("all_rd150_DupOnly.vcf",'r')
    vcf_out = pysam.VariantFile("all_dup_genotyped_12_16.vcf",'w', header=vcf_in.header)
    inside_threshold = 450 #### threshold for reads inside the signal breakpoints for that SV type
    outside_threshold = 450 #### threshold for reads outside the signal breakpoints for that SV type
    sam_iter_buffer = 600
    for var in vcf_in.fetch():
        start_pos = var.start
        end_pos = var.stop
        chrom = var.chrom

        #### breakpoint threshold dynamic on SV size
        if (end_pos - start_pos) < 300:
            inside_threshold = 150
            outside_threshold = 150
        elif 300 <= (end_pos - start_pos) < 1000:
            inside_threshold = 300
            outside_threshold = 300
        else:
            inside_threshold = 450
            outside_threshold = 450

        with open("bams.txt", "r") as bam_list:
            for f in bam_list:
                paired_support_name_list = []
                paired_lowMQ_name_list = []
                paired_support_counter = 0
                split_support_counter = 0
                bam_file = pysam.Samfile(f.strip(), "rb")
                sample = str(os.path.basename(bam_file.filename.decode('utf-8')))
                sample = sample[:sample.index("_sorted_dedup")]
                chrom_len = bam_file.get_reference_length(chrom)

                ### Check to make sure we don't have duplicated reads if svlen < sam_iter_buffer*2, min/max to make sure we aren't starting past chromosome ends
                if end_pos - start_pos <= 2*sam_iter_buffer:
                    sam_iter = bam_file.fetch(chrom, max(0, start_pos-sam_iter_buffer), min(chrom_len, end_pos+sam_iter_buffer))
                else:
                    sam_iter_start = bam_file.fetch(chrom, max(0, start_pos-sam_iter_buffer), start_pos+sam_iter_buffer)
                    sam_iter_end = bam_file.fetch(chrom, end_pos-sam_iter_buffer, min(chrom_len, end_pos+sam_iter_buffer))
                    sam_iter = chain(sam_iter_start, sam_iter_end)

                for read in sam_iter:
                    
                    #####                      #####
                    # Check split reads for signal #
                    #####                      #####

                    if (
                        read.has_tag("SA")
                        and not read.is_secondary
                        and not  read.is_qcfail
                        and not  read.is_unmapped
                        and not  read.is_duplicate
                        and not  read.is_supplementary
                        and read.mapping_quality > 0
                    ):

                        split_type = ""
                        ##### Gather info on the primary, non-supplemental part of the read
                        primary_reference = read.reference_name
                        primary_strand = not read.is_reverse  ## forward strand = True, reverse strand = False
                        primary_q_start, primary_q_end = calc_query_pos_from_cigar(read.cigarstring, primary_strand)
                        primary_start_pos = read.reference_start
                        primary_end_pos = read.reference_end

                        ##### Gather info on the supplemental part of the read attached to the SA tag
                        supp_reference = read.get_tag("SA").split(",")[0]
                        if read.get_tag("SA").split(",")[2] == "+":
                            supp_strand = True
                        else:
                            supp_strand = False
                        supp_start_pos = int(read.get_tag("SA").split(",")[1])
                        supp_q_start, supp_q_end = calc_query_pos_from_cigar(read.get_tag("SA").split(",")[3], supp_strand)
                        supp_end_pos = supp_start_pos + (supp_q_end - supp_q_start)
                        
                        ##### Create primary and supplemental mixed type lists and get split_type
                        primary_read = [primary_reference, primary_strand, primary_q_start, primary_start_pos, primary_end_pos, bam_file.get_tid(primary_reference)]
                        supp_read = [supp_reference, supp_strand, supp_q_start, supp_start_pos, supp_end_pos, bam_file.get_tid(supp_reference)]
                        
                        split_type = get_split_event_type(primary_read, supp_read)
                        
                        ##### Check if split is same as SV type and breakpoints within certain threshold of variant breakpoints

                        if var.info["SVTYPE"] == split_type:
                            left_bp, right_bp = get_breakpoints(primary_read, supp_read, split_type)
                            if left_bp in range(start_pos-outside_threshold,start_pos+inside_threshold) and right_bp in range(end_pos-inside_threshold,end_pos+outside_threshold):
                                split_support_counter += 1

                    #####                       #####
                    # Check paired reads for signal #
                    #####                       #####

                    # Here we build a list of all supporting paired reads, check initially that if one read is supported 
                    # based on paired position data w/ quality checks for that read, that if the other read fails quality
                    # checks as we iterate through reads, that pair is removed from the support list


                    if (read.query_name in paired_support_name_list):
                        if (
                            (read.mapping_quality < 1 and not read.is_supplementary)
                            or read.is_secondary
                            or read.is_qcfail
                            or read.is_unmapped
                            or read.is_duplicate
                        ):
                            paired_support_name_list.remove(read.query_name)

                    elif not (read.query_name in paired_lowMQ_name_list):
                        if(
                            read.is_qcfail
                            or read.is_unmapped
                            or read.is_duplicate
                            or read.is_secondary
                            or not read.is_paired
                            or read.is_supplementary
                        ):
                            continue

                        #### check duplication ####
                        if read.mapping_quality > 0: 
                            paired_read_q_start, paired_read_q_end = calc_query_pos_from_cigar(read.cigarstring, not read.is_reverse)
                        else:
                            paired_read_q_start = 0
                            paired_read_q_end = 0
                        paired_mate_q_start, paired_mate_q_end = calc_query_pos_from_cigar(robust_get_MC(read), not read.mate_is_reverse)
                        paired_read_end = read.reference_end
                        paired_mate_end = read.next_reference_start + (paired_mate_q_end - paired_mate_q_start)
                        
                        ##### Check duplication signal if SV = DUP
                        if var.info["SVTYPE"] == "DUP" and read.reference_id == read.next_reference_id:
                            if read.is_reverse and paired_read_end in range(start_pos-outside_threshold,start_pos+inside_threshold) and not read.mate_is_reverse:
                                if read.next_reference_start > read.reference_start and read.next_reference_start in range(end_pos-inside_threshold,end_pos+outside_threshold):
                                    if read.mapping_quality > 0:
                                        paired_support_name_list.append(read.query_name)
                                    else:
                                        paired_lowMQ_name_list.append(read.query_name)

                            if not read.is_reverse and read.reference_start in range(end_pos-inside_threshold,end_pos+outside_threshold) and read.mate_is_reverse:
                                if read.next_reference_start < read.reference_start and paired_mate_end in range(start_pos-outside_threshold,start_pos+inside_threshold):
                                    if read.mapping_quality > 0:
                                        paired_support_name_list.append(read.query_name)
                                    else:
                                        paired_lowMQ_name_list.append(read.query_name)

                paired_support_counter = len(paired_support_name_list)

                ##### Update genotypes #####
                if paired_support_counter > 0 or split_support_counter > 0:
                    var.samples[sample]['GT'] = (1,1)
                else:
                    var.samples[sample]['GT'] = (0,0)   

        #### Write variant genotypes to output SV ####                 
        vcf_out.write(var)
    vcf_out.close()

if __name__ == "__main__":
    main()