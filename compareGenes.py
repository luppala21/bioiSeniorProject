
class compareGenes:

    def main(self):

        microbiome_file_name = "RESULTS-PDAC.csv"
        microbiome_file = open(microbiome_file_name, "r")

        tissue_file_name = "RESULTS-IBS.csv"
        tissue_file = open(tissue_file_name, "r")

        # ***BE SURE to change the output file name!***
        output_file_name = "RESULTS-diseases-all-v2.csv"
        output_file = open(output_file_name, "w")
        output_file.write("gene_name,pvalue,IBS_micro_logFC,PDAC_micro_logFC\n")

        microbiome_file_list = microbiome_file.readlines()[1:]
        tissue_file_list = tissue_file.readlines()[1:]

        
        for line in tissue_file_list:

            tissue_word = line.split(",")
            print(tissue_word)
            
            tissue_logFC = float(tissue_word[3].strip())
            tissue_gene = tissue_word[0].strip()
            tissue_pvalue = (float(tissue_word[1]) 
                + float(tissue_word[2]))/2

            for item in microbiome_file_list:

                microbiome_word = item.split(",")
                #print(microbiome_word)
                
                microbiome_pvalue = (float(microbiome_word[1]) 
                    + float(microbiome_word[2]))/2
                microbiome_logFC = float(microbiome_word[3].strip())
                microbiome_gene = microbiome_word[0].strip()
                #gene_function = microbiome_word[1].strip()

                if microbiome_gene == tissue_gene:
                    pvalue = str((microbiome_pvalue + tissue_pvalue)/2)

                    output_file.write(microbiome_gene + "," 
                                    + str(microbiome_pvalue) + ","
                                    + str(tissue_logFC) + ","
                                    + str(microbiome_logFC) + "\n")
                    break

        tissue_file.close()
        microbiome_file.close()
        output_file.close()

# makes an object of the mapProbes() class to test
# the main() method
objectTest = compareGenes()
objectTest.main()