

class mapProbes:

    def main(self):

        map_file_name = "HTA_2_0.chip"
        map_file = open(map_file_name, "r")

        probe_file_name = "GSE63379.top.table-IBS.csv"
        probe_file = open(probe_file_name, "r")

        output_file_name = "RESULTS-IBS.csv"
        output_file = open(output_file_name, "w")
        output_file.write("Probe,Gene,Chr,Start,End")

        map_line_list = map_file.readlines()[1:]
        probe_file_list = probe_file.readlines()[1:]

        for item in probe_file_list:

            probe_word = item.split(",")

            probe_id = probe_word[0].strip()
            probe_chr = probe_word[6].strip()
            probe_start = probe_word[7]
            probe_end = probe_word[8]

            for line in map_line_list:

                probe_word = line.split("\t")
                probeList_id = probe_word[0]
                probe_gene = probe_word[1]

                if probe_id == probeList_id:
                    
                    output_file.write(probe_id + "," +
                                    probe_gene + "," +
                                    probe_chr + "," +
                                    probe_start + "," +
                                    probe_end)


        map_file.close()
        probe_file.close()
        output_file.close()

# makes an object of the mapProbes() class to test
# the main() method
objectTest = mapProbes()
objectTest.main()