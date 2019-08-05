#include "commander.h"
#include "functions.h"
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <stdexcept>
#include <stdint.h>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
int main(int argc, char* argv[])
{
    Commander commander;
    try
    {
        if (!commander.initialize(argc, argv)) { return -1; }
    }
    catch (const std::runtime_error& er)
    {
        std::cerr << er.what() << std::endl;
        return -1;
    }
    std::cerr << "Start processing GTF file" << std::endl;
    std::unordered_map<std::string, std::unordered_set<std::string>>
        gene_id_map = gen_gene_id_map(commander.gtf());
    std::cerr << "Start processing MSigDB file" << std::endl;
    std::vector<std::string> set_names;
    std::unordered_map<std::string, std::vector<size_t>> gene_membership =
        gen_gene_membership(commander.msigdb(), gene_id_map, set_names);
    // we don't need the gene id map after processing the msigdb
    gene_id_map.clear();
    std::cerr << "Start processing bim and summary statistic file" << std::endl;
    std::unordered_map<std::string, std::string> snp_id_map = get_snps(
        commander.target(), commander.sumstat(), commander.sumstat_snp());
    if (snp_id_map.empty() || snp_id_map.size() == 0)
    {
        std::cerr << "Error: No SNPs left!" << std::endl;
        return -1;
    }
    std::cerr << "Obtain eQTL p-values" << std::endl;
    std::unordered_map<std::string, std::vector<size_t>> snp_p =
        gen_binary_pathway_member(
            commander.eqtl(), commander.eqtl_snp(), commander.eqtl_gene(),
            commander.eqtl_pvalue(), snp_id_map, gene_membership,
            commander.threshold(), set_names.size());
    generate_snp_sets(snp_p, set_names, commander.threshold(), commander.out());
    std::cerr << "Completed" << std::endl;
    return 0;
}
