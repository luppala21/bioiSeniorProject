

class mapProbes:

    def main(self):

        map_file_name = "HTA_2_0.chip"
        map_file = open(map_file_name, "r")

        probe_file_name = "GSE63379.top.table-IBS_tissue.tsv"
        probe_file = open(probe_file_name, "r")

        output_file_name = "RESULTS-IBS.csv"
        output_file = open(output_file_name, "w")
        output_file.write("ID,adj.P.val,P.Value,t,B,logFC,Gene.Symbol,Gene.Function\n")

        map_line_list = map_file.readlines()[1:]
        probe_file_list = probe_file.readlines()[1:]

        for item in probe_file_list:

            probe_word = item.split("\t")
            #print(probe_word)

            probe_id = probe_word[0].strip()
            probe_chr = probe_word[6].strip()

            for line in map_line_list:

                map_word = line.split("\t")
                #print(map_word)
                
                mapList_id = map_word[0]
                map_gene = map_word[1].strip()
                map_geneFunction = map_word[2]

                if probe_id == mapList_id:
                    
                    output_file.write(probe_id + "\t" +
                                probe_word[1] + "\t" +
                                probe_word[2] + "\t" +
                                probe_word[3] + "\t" +
                                probe_word[4] + "\t" +
                                map_gene + "\t" +
                                map_geneFunction + "\n")


        map_file.close()
        probe_file.close()
        output_file.close()

# makes an object of the mapProbes() class to test
# the main() method
objectTest = mapProbes()
objectTest.main()