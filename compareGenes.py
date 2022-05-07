class compareGenes:

    def main(self):

        # name of the microbiome gene list input file
        microbiome_file_name = "RESULTS-PDAC.csv"
        microbiome_file = open(microbiome_file_name, "r")

        # name of the tissue gene list input file
        tissue_file_name = "RESULTS-IBS.csv"
        tissue_file = open(tissue_file_name, "r")

        # name of output file which will contain
        # genes common to both input files
        # termed the "intersection" of the files
        output_file_name = "RESULTS-diseases-all-v2.csv"
        output_file = open(output_file_name, "w")

        # header of output file
        output_file.write("gene_name,pvalue,IBS_micro_logFC,PDAC_micro_logFC\n")

        microbiome_file_list = microbiome_file.readlines()[1:]
        tissue_file_list = tissue_file.readlines()[1:]

        
        # iterates through the lists and outputs a list
        # which contains genes which are common to both
        for line in tissue_file_list:

            tissue_word = line.split(",")
            print(tissue_word)
            
            tissue_logFC = float(tissue_word[3].strip())
            tissue_gene = tissue_word[0].strip()
            #tissue_pvalue = (float(tissue_word[1]) 
            #    + float(tissue_word[2]))/2

            for item in microbiome_file_list:

                microbiome_word = item.split(",")
                #print(microbiome_word)
                
                microbiome_pvalue = (float(microbiome_word[1]) 
                    + float(microbiome_word[2]))/2
                microbiome_logFC = float(microbiome_word[3].strip())
                microbiome_gene = microbiome_word[0].strip()
                #gene_function = microbiome_word[1].strip()

                if microbiome_gene == tissue_gene:
                    #pvalue = str((microbiome_pvalue + tissue_pvalue)/2)

                    output_file.write(microbiome_gene + "," 
                                    + str(microbiome_pvalue) + ","
                                    + str(tissue_logFC) + ","
                                    + str(microbiome_logFC) + "\n")
                    break

        tissue_file.close()
        microbiome_file.close()
        output_file.close()

# makes an object of the compareGenes() class to test
# the main() method
objectTest = compareGenes()
objectTest.main()