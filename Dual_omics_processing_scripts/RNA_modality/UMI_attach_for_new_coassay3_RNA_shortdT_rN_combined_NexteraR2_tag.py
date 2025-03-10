import subprocess
import sys
from Levenshtein import distance
import gzip
import scipy.spatial.distance
from multiprocessing import Pool
from functools import partial
import pickle
import copy


def UMI_attach_read2_barcode_list(sample, input_folder, output_folder, ligation_1st_barcode_list, ligation_2nd_barcode_list, mismatch_rate = 1):

    # sample = 'sciATAC3_EXP10_01'
    Read1 = input_folder + "/" + sample + ".R1.fastq.gz"
    Read2 = input_folder + "/" + sample + ".R3.fastq.gz"
    Read3 = input_folder + "/" + sample + ".R2.fastq.gz"
    
    output_file1 = output_folder + "/" + sample + "_dT.R2.fastq.gz"
    output_file2 = output_folder + "/" + sample + "_rN.R2.fastq.gz"

    f1 = gzip.open(Read1)
    f2 = gzip.open(Read2)
    I5_reads = gzip.open(Read3)
    
    f3 = gzip.open(output_file1, 'wb')
    f4 = gzip.open(output_file2, 'wb')

    line1 = f1.readline()
    line2 = f2.readline()
    line3 = I5_reads.readline()
    
    total_PE_line = 0
    dT_RNA_PE_line = 0
    rN_RNA_PE_line = 0
    ligation_matched_PE_line  = 0
    
    ###generate the error list of polyT and ME sequence
    ME = "GAGATGTGTA"
    ME_list = list(ME)
    ME_error_list = set()
    for first_error in range(len(ME)):
        for each_base in ["A", "G", "C", "T", "N"]:
            ME_list_first_error = copy.deepcopy(ME_list)
            ME_list_first_error[first_error] = each_base
            ME_list_first_error_seq = "".join(ME_list_first_error)
            ME_error_list.add(ME_list_first_error_seq)
            if first_error+1 < len(ME):
                for second_error in range(first_error+1,len(ME)):
                    for each_base_again in ["A", "G", "C", "T", "N"]:
                        ME_list_second_error = copy.deepcopy(ME_list_first_error)
                        ME_list_second_error[second_error] = each_base_again
                        ME_list_second_error_seq = "".join(ME_list_second_error)
                        ME_error_list.add(ME_list_second_error_seq)                    
    ME_error_list = list(ME_error_list)
    
    polyT = "TTTTTTTTTT"
    polyT_list = list(polyT)
    polyT_error_list = set()
    for first_error in range(len(polyT)):
        for each_base in ["A", "G", "C", "T", "N"]:
            polyT_list_first_error = copy.deepcopy(polyT_list)
            polyT_list_first_error[first_error] = each_base
            polyT_list_first_error_seq = "".join(polyT_list_first_error)
            polyT_error_list.add(polyT_list_first_error_seq)               
    polyT_error_list = list(polyT_error_list)
    
    cocapture_constant1 = "CAAGTTGATA"
    constant_mismatch_list1 = set()
    for each_pos in range(len(cocapture_constant1)):
        ori_list = list(cocapture_constant1)
        for each_letter in ["A", "G", "C", "T", "N"]:
            ori_list[each_pos] = each_letter
            constant_mismatch_list1.add("".join(ori_list))
    constant_mismatch_list1 = list(constant_mismatch_list1)

    while (line1):
        total_PE_line += 1
        line1_header = line1
        line2_header = line2
        line1 = f1.readline()
        line2 = f2.readline()
        line3 = I5_reads.readline()
            #print("read1: ", line1)
            # first check if the ligation barcode match with the barcode
        tmp_lig1 = line1[0:10].decode()
        tmp_lig2 = line3[0:10].decode()
        RNA_primer_identifier = line1[40:42].decode()  ###should be TG
        ban_seq = line1[50:56].decode() ###should be randomN sequence if is from rN
        potential_ME_seq = line1[40:50].decode()

        ## check if the ligation is successful
        if tmp_lig2 in ligation_2nd_barcode_list and tmp_lig1 in ligation_1st_barcode_list:
            ligation_bc2_match = ligation_2nd_barcode_list[tmp_lig2]
            ligation_bc1_match = ligation_1st_barcode_list[tmp_lig1]
            ligation_matched_PE_line += 1
            
            # get the potential dT/cocapture pattern
            polyT_seq = line1[50:60].decode()
            
            ## check if the read is from dT priming
            if polyT_seq in polyT_error_list and RNA_primer_identifier == "TG" and potential_ME_seq not in ME_error_list:
                dT_RNA_PE_line += 1
                UMI = line1[42:50].decode()
                first_line = '@' + ligation_bc1_match + ligation_bc2_match + ',' + UMI + ',' +line2_header[1:].decode()
                f3.write(first_line.encode())
                
                second_line = line2
                f3.write(second_line)

                line1_third = f1.readline()
                third_line = f2.readline()
                f3.write(third_line)

                line1_fourth = f1.readline()
                four_line = f2.readline()
                f3.write(four_line)

                line1 = f1.readline()
                line2 = f2.readline()
            
            ###check if is from randomN priming
            elif RNA_primer_identifier == "TG" and polyT_seq not in polyT_error_list and \
            polyT_seq not in constant_mismatch_list1 \
            and ban_seq != "TTTTTT" and potential_ME_seq not in ME_error_list:
                rN_RNA_PE_line += 1
                UMI = line1[42:50].decode()
                first_line = '@' + ligation_bc1_match + ligation_bc2_match + ',' + UMI + ',' +line2_header[1:].decode()
                
                f4.write(first_line.encode())
                
                second_line = line2
                f4.write(second_line)

                line1_third = f1.readline()
                third_line = f2.readline()
                f4.write(third_line)

                line1_fourth = f1.readline()
                four_line = f2.readline()
                f4.write(four_line)

                line1 = f1.readline()
                line2 = f2.readline()
                
            else:
                line1 = f1.readline()
                line1 = f1.readline()
                line1 = f1.readline()
                line2 = f2.readline()
                line2 = f2.readline()
                line2 = f2.readline()
        else:      
            line1 = f1.readline()
            line1 = f1.readline()
            line1 = f1.readline()
            line2 = f2.readline()
            line2 = f2.readline()
            line2 = f2.readline()
        line3 = I5_reads.readline() 
        line3 = I5_reads.readline()
        line3 = I5_reads.readline()
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    I5_reads.close()

    print("sample name: %s, total PE line: %f, ligation barcodes matched PE line: %f, \
          PE line from dT primer: %f, dT filter rate: %f, \
          PE line from rN primer: %f, rN filter rate: %f" 
           % (sample, total_PE_line, ligation_matched_PE_line, \
              dT_RNA_PE_line, float(dT_RNA_PE_line)/float(total_PE_line),\
              rN_RNA_PE_line, float(rN_RNA_PE_line)/float(total_PE_line)))


def attach_UMI_files(input_folder, sampleID, output_folder, ligation_barcode1_file, ligation_barcode2_file, core):
    
    init_message = '''
    --------------------------start attaching UMI-----------------------------
    input folder: %s
    sample ID: %s
    output_folder: %s
    ligation barcode1 file: %s
    ligation barcode2 file: %s
    ___________________________________________________________________________
    ''' %(input_folder, sampleID, output_folder, ligation_barcode1_file, ligation_barcode2_file)
    
    print(init_message)
    
    
    print("Load ligation barcodes dictionary...")

    barcodes = open(ligation_barcode1_file, "rb")
    ligation_1st_barcode_list = pickle.load(barcodes)
    barcodes.close()

    barcodes = open(ligation_barcode2_file, "rb")
    ligation_2nd_barcode_list = pickle.load(barcodes)
    barcodes.close()

    #for each sample in the sample list, use the read1 file, read2 file, output file
    # and barcode_list to run UMI_attach_read2_barcode_list
    sample_file = open(sampleID)
    sample_list = []
    for line in sample_file:
        sample = line.strip()
        sample_list.append(sample)
    sample_file.close()
    
    # parallele for the functions
    p = Pool(processes = int(core))
    #print("Processing core number: ", core_number)
    func = partial(UMI_attach_read2_barcode_list, input_folder = input_folder, output_folder=output_folder, ligation_1st_barcode_list=ligation_1st_barcode_list, ligation_2nd_barcode_list=ligation_2nd_barcode_list, mismatch_rate = 1)
    #sciRNAseq_count(sample, input_folder, exons, genes, gene_end)
    result = p.map(func, sample_list)
    p.close()
    p.join()
    
    #print the completion message
    com_message = '''~~~~~~~~~~~~~~~UMI attachment done~~~~~~~~~~~~~~~~~~'''
    print(com_message)
    
if __name__ == "__main__":
    input_folder = sys.argv[1]
    sampleID = sys.argv[2]
    output_folder = sys.argv[3]
    ligation_barcode1_file = sys.argv[4]
    ligation_barcode2_file = sys.argv[5]
    core=sys.argv[6]
    attach_UMI_files(input_folder, sampleID, output_folder, ligation_barcode1_file, ligation_barcode2_file, core)
