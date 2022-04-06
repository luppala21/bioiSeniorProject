

class compareGenes:

    def main(self):

        microbiome_file_name = "GSE14841.top.table-IBS_microbiome.tsv"
        microbiome_file = open(microbiome_file_name, "r")

        tissue_file_name = "RESULTS-IBS.tsv"
        tissue_file = open(tissue_file_name, "r")

        output_file_name = "RESULTS-IBS-geneCompare.csv"
        output_file = open(output_file_name, "w")
        output_file.write("gene_name,gene_function")

        microbiome_file_list = microbiome_file.readlines()[1:]
        tissue_file_list = tissue_file.readlines()[1:]

        for item in microbiome_file_list:

            microbiome_word = item.split("\t")
            #print(microbiome_word)
            
            microbiome_gene = microbiome_word[6].strip()
            gene_function = microbiome_word[7].strip()

            for line in tissue_file_list:

                tissue_word = line.split("\t")
                #print(tissue_word)
                
                tissue_gene = tissue_word[5].strip()

                if microbiome_gene == tissue_gene:

                    output_file.write(microbiome_gene + "," + gene_function + "\n")

            
            


        tissue_file.close()
        microbiome_file.close()
        output_file.close()

# makes an object of the mapProbes() class to test
# the main() method
objectTest = compareGenes()
objectTest.main()