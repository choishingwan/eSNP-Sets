#ifndef COMMANDER_H
#define COMMANDER_H

#include "misc.hpp"
#include <getopt.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

class Commander
{
public:
    Commander();
    virtual ~Commander();
    bool initialize(int argc, char* argv[]);
    std::vector<double> threshold() const { return m_p_threshold; }
    std::string eqtl() const { return m_eqtl_name; }
    std::string sumstat() const { return m_sumstat_name; }
    std::string msigdb() const { return m_msigdb_name; }
    std::string gtf() const { return m_gtf_name; }
    std::string target() const { return m_target_name; }
    std::string eqtl_snp() const { return m_e_snp_id; }
    std::string eqtl_gene() const { return m_e_gene_id; }
    std::string eqtl_pvalue() const { return m_e_pvalue; }
    std::string sumstat_snp() const { return m_s_snp_id; }
    std::string out() const { return m_out; }

private:
    std::vector<double> m_p_threshold;
    std::string m_eqtl_name;
    std::string m_sumstat_name;
    std::string m_msigdb_name;
    std::string m_gtf_name;
    std::string m_target_name;
    std::string m_e_snp_id = "variant_id";
    std::string m_e_gene_id = "gene_id";
    std::string m_e_pvalue = "";
    std::string m_s_snp_id = "";
    std::string m_out = "out";
    std::string m_version = "0.0.1";
    std::string m_date = "2019-08-02";
    void usage();
};

#endif // COMMANDER_H
