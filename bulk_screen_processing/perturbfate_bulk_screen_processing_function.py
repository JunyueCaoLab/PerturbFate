import pandas as pd
import gzip
from Bio.Seq import Seq

def trisci_bulk_screen_processing(fq_R1_dir, fq_R2_dir, fq_I5_dir, i5_seq_for_demultiplexing, seq_gene_table_dir, output_count_csv_dir):    
    
    ###generate the white list
    screen_ref_table = pd.read_csv(seq_gene_table_dir, header=0)
    
    ###add the column of gRNA combination
    screen_ref_table["Barcode_Sequence"] = screen_ref_table["sgRNA_1_seq"] + screen_ref_table["sgRNA_2_seq"]
    
    ###build the sgRNA whitelist
    print("Building the sgRNA whitelist...")
    seq_whitelist_dict = {}

    for each_gene in screen_ref_table.index:
        gRNA1 = str(screen_ref_table.loc[each_gene, "sgRNA_1_seq"])
        gRNA2 = str(screen_ref_table.loc[each_gene, "sgRNA_2_seq"])

        gRNA1_mismatch_list = set()
        gRNA2_mismatch_list = set()

        for each_pos in range(len(gRNA1)):
            ori_gRNA1_list = list(gRNA1)
            for each_letter in ["A", "G", "C", "T", "N"]:
                ori_gRNA1_list[each_pos] = each_letter
                gRNA1_mismatch_list.add("".join(ori_gRNA1_list))
        gRNA1_mismatch_list = list(gRNA1_mismatch_list)

        for each_pos in range(len(gRNA2)):
            ori_gRNA2_list = list(gRNA2)
            for each_letter in ["A", "G", "C", "T", "N"]:
                ori_gRNA2_list[each_pos] = each_letter
                gRNA2_mismatch_list.add("".join(ori_gRNA2_list))
        gRNA2_mismatch_list = list(gRNA2_mismatch_list)

        gRNA_combination = []
        for each_gRNA1 in gRNA1_mismatch_list:
            for each_gRNA2 in gRNA2_mismatch_list:
                gRNA_combination.append(each_gRNA1 + each_gRNA2)

        seq_whitelist_dict.update(dict(zip(gRNA_combination, [gRNA1+gRNA2]*len(gRNA_combination))))
    
    ###start processing reads
    print("Start processing reads...")
    fq_R1_dir=fq_R1_dir
    fq_R2_dir=fq_R2_dir
    fq_I5_dir=fq_I5_dir
    
    gene_seq_dict = dict(zip(screen_ref_table["gene_symbol"], screen_ref_table["Barcode_Sequence"]))
    seq_count_dict = dict(zip(screen_ref_table["Barcode_Sequence"], [0]*len(screen_ref_table["Barcode_Sequence"])))

    fq_R1 = gzip.open(fq_R1_dir, "rb")
    fq_R2 = gzip.open(fq_R2_dir, "rb")
    fq_I5 = gzip.open(fq_I5_dir, "rb")

    line1 = fq_R1.readline()
    line2 = fq_R2.readline()
    line3 = fq_I5.readline()

    ori_R1_constant = "AAGGACGAAACACCG"                      
    ori_R2_constant = "GCCACTTTTTCAAGT"

    constant_R1_mismatch_list = set()                      
    constant_R2_mismatch_list = set()

    for each_pos in range(len(ori_R1_constant)):
        ori_R1_list = list(ori_R1_constant)
        for each_letter in ["A", "G", "C", "T", "N"]:
            ori_R1_list[each_pos] = each_letter
            constant_R1_mismatch_list.add("".join(ori_R1_list))
    constant_R1_mismatch_list = list(constant_R1_mismatch_list)

    for each_pos in range(len(ori_R2_constant)):
        ori_R2_list = list(ori_R2_constant)
        for each_letter in ["A", "G", "C", "T", "N"]:
            ori_R2_list[each_pos] = each_letter
            constant_R2_mismatch_list.add("".join(ori_R2_list))
    constant_R2_mismatch_list = list(constant_R2_mismatch_list)

    total_reads = 0
    single_sample_total_reads = 0
    constant_region_matched_reads = 0                      
    correct_gRNA_reads = 0                      

    while line1:
        total_reads += 1                  
        line1 = fq_R1.readline()
        line2 = fq_R2.readline()
        line3 = fq_I5.readline()
        constant_R1 = line1.decode()[10:25]
        constant_R2 = line2.decode()[10:25]
        i5 = line3.decode().strip()
        
        if total_reads%1000000==0:
            print(str(total_reads) + " reads have been processed.")
        
        if i5 == i5_seq_for_demultiplexing:
            single_sample_total_reads += 1

            if constant_R1 in constant_R1_mismatch_list and constant_R2 in constant_R2_mismatch_list:                  
                constant_region_matched_reads += 1
                gRNA_seq1 = line1.decode()[24:44]
                gRNA_seq2_rc = line2.decode()[83:103]

                ###perform rc on sgRNA seq2
                gRNA_seq2_rc_seq = Seq(gRNA_seq2_rc)
                gRNA_seq2 = str(gRNA_seq2_rc_seq.reverse_complement())


                #allow 1 mismatch of each sgRNA
                if (gRNA_seq1 + gRNA_seq2) in seq_whitelist_dict.keys():

                    correct_gRNA_reads += 1
                    corrected_gRNA_combinations = seq_whitelist_dict[(gRNA_seq1 + gRNA_seq2)]
                    seq_count_dict[corrected_gRNA_combinations] += 1


        line1 = fq_R1.readline()
        line2 = fq_R2.readline()
        line3 = fq_I5.readline()
        line1 = fq_R1.readline()
        line2 = fq_R2.readline()
        line3 = fq_I5.readline()
        line1 = fq_R1.readline()
        line2 = fq_R2.readline()
        line3 = fq_I5.readline()

    print("total reads processed: %s\ntotal reads of this sample: %s\nconstant region matched reads processed: %s\ngRNA reads processed: %s" % \
          (total_reads, single_sample_total_reads, constant_region_matched_reads, correct_gRNA_reads))                      

    fq_R1.close()
    fq_R2.close()
    fq_I5.close()
    
    print("formatting the result...")
    
    #for gRNA_seq, counts in seq_count_dict.items():
    count_df = pd.DataFrame.from_dict(seq_count_dict, orient="index", columns=["Lib_counts"])
    count_df["Barcode_Sequence"] = count_df.index
    screen_ref_table.index = screen_ref_table["Barcode_Sequence"]
    count_output_table = pd.concat([screen_ref_table, count_df], axis=1, join="inner")
    count_output_table.to_csv(output_count_csv_dir, header = True, index = False)
    
    print("All done!")
    
    return count_output_table