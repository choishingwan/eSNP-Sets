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
    std::unordered_map<std::string, std::unordered_set<std::string>>
        gene_id_map = gen_gene_id_map(commander.gtf());
    std::vector<std::string> set_names;
    std::unordered_map<std::string, std::vector<uint64_t>> gene_membership =
        gen_gene_membership(commander.msigdb(), gene_id_map, set_names);
    // we don't need the gene id map after processing the msigdb
    gene_id_map.clear();
    std::unordered_map<std::string, std::string> snp_id_map = get_snps(
        commander.target(), commander.sumstat(), commander.sumstat_snp());
    std::vector<std::vector<double>> snp_p = gen_pathway_member(
        commander.eqtl(), commander.eqtl_snp(), commander.eqtl_gene(),
        commander.eqtl_pvalue(), snp_id_map, gene_membership, set_names.size());
    return 0;
}
