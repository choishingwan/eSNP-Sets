#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "gzstream.h"
#include "misc.hpp"
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>


std::unordered_map<std::string, std::unordered_set<std::string>>
gen_gene_id_map(const std::string& gtf_name)
{
    std::unordered_map<std::string, std::unordered_set<std::string>> result;
    bool is_gz = gtf_name.substr(gtf_name.find_last_of(".") + 1) == ("gz");
    std::ifstream gtf;
    GZSTREAM_NAMESPACE::igzstream gz_gtf;
    std::string error_message;
    if (is_gz)
    {
        gz_gtf.open(gtf_name.c_str());
        if (!gz_gtf.good())
        {
            error_message = "Error: Cannot open gziped GTF file - " + gtf_name;
            throw std::runtime_error(error_message);
        }
    }
    else
    {
        gtf.open(gtf_name.c_str());
        if (!gtf.is_open())
        {
            error_message = "Error: Cannot open GTF file - " + gtf_name;
            throw std::runtime_error(error_message);
        }
    }

    std::string line;
    std::vector<std::string> token, attribute, extract;
    bool has_gene_id = false, has_gene_name = false;
    std::string gene_id, gene_name;
    while ((is_gz && std::getline(gz_gtf, line))
           || (!is_gz && std::getline(gtf, line)))
    {
        misc::trim(line);
        if (line.empty()) continue;
        if (line.front() == '#') continue;
        token = misc::split(line, "\t");
        // now get the attributes
        if (token.size() != 9)
        {
            error_message =
                "Error: Malformed GTF file. GTF file should contain 9 columns: "
                + line;
            throw std::runtime_error(error_message);
        }
        attribute = misc::split(token[8], ";");
        has_gene_id = false;
        has_gene_name = false;
        for (auto&& a : attribute)
        {
            extract = misc::split(a, " ");
            if (extract.size() != 2)
            {
                error_message = "Error: Malformed attribute in GTF: " + line;
                throw std::runtime_error(error_message);
            }
            if (extract.front() == "gene_id")
            {
                has_gene_id = true;
                gene_id = extract.back();
            }
            else if (extract.front() == "gene_name")
            {
                has_gene_name = true;
                gene_name = extract.back();
            }
            if (has_gene_id && has_gene_name) break;
        }
        if (has_gene_id && has_gene_name)
        {
            // allow 1 to many mapping
            result[gene_name].insert(gene_id);
        }
    }
    if (is_gz) { gz_gtf.close(); }
    else
    {
        gtf.close();
    }
    std::cerr << result.size() << " gene identified in the GTF file"
              << std::endl;
    return result;
}

std::unordered_map<std::string, std::vector<uint64_t>> gen_gene_membership(
    const std::string& msigdb_name,
    const std::unordered_map<std::string, std::unordered_set<std::string>>&
        gene_map,
    std::vector<std::string>& set_name)
{
    std::unordered_map<std::string, std::vector<uint64_t>> result;
    std::ifstream msigdb;
    msigdb.open(msigdb_name);
    std::string error_message = "";
    if (!msigdb.is_open())
    {
        error_message = "Error: Cannot open MSigDB file - " + msigdb_name;
        throw std::runtime_error(error_message);
    }
    std::string line;
    std::vector<std::string> token;
    size_t num_sets = 0;
    while (std::getline(msigdb, line))
    {
        misc::trim(line);
        num_sets += !line.empty();
    }
    msigdb.clear();
    msigdb.seekg(0, std::ios::beg);
    // size of the uint64_t vector
    // we + 1 to the nunber of sets as we want to include a set that include all
    // eQTL signal
    size_t flag_size = ((num_sets + 1) / 64) + 1;
    // idx starts at 1, because 0 is reserved for background set
    size_t idx = 1;
    set_name.push_back("Background");
    while (std::getline(msigdb, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        token = misc::split(line);
        if (token.size() < 2)
        {
            error_message = "Error: Malformed MSigDB file. It should contain "
                            "at least 2 columns: "
                            + line;
            throw std::runtime_error(error_message);
        }
        if (result.find(token[0]) != result.end())
        {
            std::cerr << "Warning: Duplicated Pathway - " << token[0]
                      << std::endl;
            std::cerr << "Only the first instance will be used" << std::endl;
            continue;
        }
        set_name.push_back(token[0]);
        for (size_t i = 1; i < token.size(); ++i)
        {
            // first, translate the gene name to gene id
            if (gene_map.find(token[i]) == gene_map.end()) continue;
            for (auto&& gene_id : gene_map.at(token[i]))
            {
                if (result.find(gene_id) != result.end())
                { result[gene_id][idx / 64] |= 1ULL << ((idx) % 64); }
                else
                {
                    std::vector<uint64_t> tmp(flag_size, 0ULL);
                    tmp[0 / 64] |= 1ULL;
                    tmp[idx / 64] |= 1ULL << ((idx) % 64);
                    result[gene_id] = tmp;
                }
            }
        }
        ++idx;
    }
    msigdb.close();
    std::cerr << num_sets << " Set(s) found in the MSigDB file" << std::endl;
    return result;
}

std::unordered_map<std::string, std::string>
get_snps(const std::string& target, const std::string& sumstat_name,
         const std::string& snp_id)
{
    std::unordered_map<std::string, std::string> result;
    // first check the target

    std::string error_message = "";
    std::ifstream sumstat;
    GZSTREAM_NAMESPACE::igzstream gz_sumstat;
    bool is_gz =
        sumstat_name.substr(sumstat_name.find_last_of(".") + 1) == ("gz");
    if (is_gz)
    {
        gz_sumstat.open(sumstat_name.c_str());
        if (!gz_sumstat.good())
        {
            error_message = "Error: Cannot open gziped summary statistic file: "
                            + sumstat_name;
            throw std::runtime_error(error_message);
        }
    }
    else
    {
        sumstat.open(sumstat_name.c_str());
        if (!sumstat.is_open())
        {
            error_message =
                "Error: Cannot open summary statistic file: " + sumstat_name;
            throw std::runtime_error(error_message);
        }
    }
    std::string line;
    if (is_gz)
        std::getline(gz_sumstat, line);
    else
        std::getline(sumstat, line);
    misc::trim(line);
    if (line.empty())
        throw std::runtime_error(
            "Error: Header of summary statistic file cannot be empty!");
    std::vector<std::string> token = misc::split(line);
    size_t snp_idx = 0;
    for (; snp_idx < token.size(); ++snp_idx)
    {
        if (token[snp_idx] == snp_id) { break; }
    }
    if (snp_idx == token.size() || token.size() == 0)
    {
        error_message = "Error: Cannot find the SNP ID column - " + snp_id
                        + "\nHeader is:\n" + line;
        throw std::runtime_error(error_message);
    }
    std::unordered_set<std::string> snp_in_sumstat;
    while ((is_gz && getline(gz_sumstat, line))
           || (!is_gz && getline(sumstat, line)))
    {
        misc::trim(line);
        if (line.empty()) continue;
        token = misc::split(line);
        if (token.size() < snp_idx)
        {
            error_message =
                "Error: Malformed summary statistic file. Should have at least "
                + std::to_string(snp_idx) + " columns";
            throw std::runtime_error(error_message);
        }
        snp_in_sumstat.insert(token[snp_idx]);
    }
    if (is_gz)
        gz_sumstat.close();
    else
        sumstat.close();
    std::cerr << snp_in_sumstat.size()
              << " SNPs found in the summary statistic file" << std::endl;
    // now process the bim file, where we have a well defined format
    std::ifstream bim;
    bool is_bim = target.substr(target.find_last_of(".") + 1) == ("bim");
    if (is_bim) { bim.open(target.c_str()); }
    else
    {
        bim.open(std::string(target + ".bim").c_str());
    }
    if (!bim.is_open())
    {
        error_message = "Error: Cannot open bim file: ";
        if (is_bim)
            error_message.append(target);
        else
            error_message.append(target + ".bim");
        throw std::runtime_error(error_message);
    }
    size_t num_not_found = 0;
    size_t num_line = 0;
    while (std::getline(bim, line))
    {
        misc::trim(line);
        if (line.empty()) continue;
        ++num_line;
        token = misc::split(line);
        if (token.size() != 6)
        {
            error_message =
                "Error: Malformed bim file, should contain 6 columns: " + line;
            throw std::runtime_error(error_message);
        }
        if (snp_in_sumstat.find(token[1]) != snp_in_sumstat.end())
        {
            // this snp is found in the summary statistic
            // we assume the chromosome ID doesn't starts with CHR
            std::string id = token[0] + "_" + token[3];
            if (result.find(id) != result.end())
            {
                std::cerr << "Warning: Duplicated SNP ID generated - " << id
                          << std::endl;
            }
            result[id] = token[1];
        }
        else
        {
            ++num_not_found;
        }
    }
    bim.close();
    // now we can read in the snp id from the summary statistic and generate the
    // location information
    std::cerr << num_line << " SNPs in bim file" << std::endl;
    std::cerr << num_not_found
              << " SNPs in bim file not found in the summary statistic file"
              << std::endl;
    std::cerr << result.size()
              << " SNPs common to the bim and summary statistic file found"
              << std::endl;

    return result;
}

std::unordered_map<std::string, std::vector<double>> gen_pathway_member(
    const std::string& eqtl_name, const std::string& eqtl_snp_id,
    const std::string& eqtl_gene_id, const std::string& eqtl_p,
    const std::unordered_map<std::string, std::string>& valid_snps,
    const std::unordered_map<std::string, std::vector<uint64_t>>&
        gene_membership,
    const size_t num_set)
{

    std::unordered_map<std::string, std::vector<double>> snps;
    bool is_gz = eqtl_name.substr(eqtl_name.find_last_of(".") + 1) == ("gz");
    std::ifstream eqtl;
    GZSTREAM_NAMESPACE::igzstream gz_eqtl;
    std::string error_message = "";
    std::string line;
    if (is_gz)
    {
        gz_eqtl.open(eqtl_name.c_str());
        if (!gz_eqtl.good())
        {
            error_message = "Error: Cannot open gziped eQTL file: " + eqtl_name;
            throw std::runtime_error(error_message);
        }
        std::getline(gz_eqtl, line);
    }
    else
    {
        eqtl.open(eqtl_name.c_str());
        if (!eqtl.is_open())
        {
            error_message = "Error: Cannot open eQTL file: " + eqtl_name;
            throw std::runtime_error(error_message);
        }
        std::getline(eqtl, line);
    }
    misc::trim(line);
    if (line.empty())
    {
        error_message = "Error: Header of eQTL file is empty";
        throw std::runtime_error(error_message);
    }

    size_t snp_idx = 0, p_idx = 0, gene_idx = 0;
    bool snp_found = false, p_found = false, gene_found = false;
    std::vector<std::string> token = misc::split(line);
    for (size_t i = 0; i < token.size(); ++i)
    {
        if (token[i] == eqtl_snp_id)
        {
            snp_idx = i;
            snp_found = true;
        }
        if (token[i] == eqtl_gene_id)
        {
            gene_idx = i;
            gene_found = true;
        }
        if (token[i] == eqtl_p)
        {
            p_idx = i;
            p_found = true;
        }
    }
    if (!snp_found || !p_found || !gene_found)
    {
        error_message =
            "Error: Required columns not found in the eQTL file: " + line;
        throw std::runtime_error(error_message);
    }
    size_t max_idx = 0;
    if (snp_idx > max_idx) max_idx = snp_idx;
    if (gene_idx > max_idx) max_idx = gene_idx;
    if (p_idx > max_idx) max_idx = p_idx;
    // now process the eqtl file
    size_t total_entry = 0, included_entry = 0;
    std::vector<std::string> snp_id;
    std::string cur_id, gene_name;
    double p_value = 2.0;
    while ((is_gz && getline(gz_eqtl, line)) || (!is_gz && getline(eqtl, line)))
    {
        misc::trim(line);
        if (line.empty()) continue;
        total_entry++;
        token = misc::split(line);
        if (max_idx >= token.size())
        {
            error_message = "Error: eQTL file does not contain enough columns!";
            throw std::runtime_error(error_message);
        }
        snp_id = misc::split(token[snp_idx], "_");
        // only use chr and bp to match
        if (snp_id.size() < 2)
        {
            error_message =
                "Error: Malform SNP ID. Should have format CHR_BP_A1_A2";
            throw std::runtime_error(error_message);
        }
        cur_id = snp_id[0] + "_" + snp_id[1];
        auto&& res = valid_snps.find(cur_id);
        if (res == valid_snps.end()) continue;
        // found
        included_entry++;
        gene_name = token[gene_idx];
        try
        {
            p_value = misc::convert<double>(token[p_idx]);
        }
        catch (const std::runtime_error&)
        {
            // check if it is NA
            std::transform(token[p_idx].begin(), token[p_idx].end(),
                           token[p_idx].begin(), ::toupper);
            if (token[p_idx] != "NAN" && token[p_idx] != "NA"
                && token[p_idx] != "NULL")
            {
                error_message =
                    "Error: Non-numeric p-value: " + misc::to_string(p_value);
                throw std::runtime_error(error_message);
            }
        }
        auto&& gene_info = gene_membership.find(gene_name);
        // Gene isn't found in MSigDB, so should be consider as part of the
        // background
        bool background_only = (gene_info == gene_membership.end());
        auto&& idx = snps.find(cur_id);
        if (idx == snps.end())
        {
            // num_set should contain the background here
            // use res->second for the RS ID instead of the eQTL ID
            snps[res->second] = std::vector<double>(num_set, 2.0);
        }

        if (background_only)
        {
            double cur_p = snps[res->second][0];
            if (cur_p > p_value) { snps[res->second][0] = p_value; }
        }
        else
        {
            // loopthrough gene_info
            for (size_t set_idx = 0; set_idx < num_set; ++set_idx)
            {
                if ((gene_info->second[(set_idx) / 64] >> ((set_idx) % 64)) & 1)
                {
                    // this is the set
                    double cur_p = snps[res->second][set_idx];
                    if (cur_p > p_value) snps[res->second][set_idx] = p_value;
                }
            }
        }
    }
    if (is_gz)
        gz_eqtl.close();
    else
        eqtl.close();
    // SNP matrix is way too big.
    // now we have SNPs, which contain the p-value of each SNP to each set
    std::cerr << snps.size() << " SNPs remaining after reading the eQTL file"
              << std::endl;
    return snps;
}


size_t calculate_column(const double& pvalue,
                        const std::vector<double>& barlevels, double& pthres)
{
    for (size_t i = 0; i < barlevels.size(); ++i)
    {
        if (pvalue < barlevels[i]
            || misc::logically_equal(pvalue, barlevels[i]))
        {
            pthres = barlevels[i];
            return i;
        }
    }
    pthres = 1.0;
    return barlevels.size();
}

std::unordered_map<std::string, std::vector<uint64_t>>
gen_binary_pathway_member(
    const std::string& eqtl_name, const std::string& eqtl_snp_id,
    const std::string& eqtl_gene_id, const std::string& eqtl_p,
    const std::unordered_map<std::string, std::string>& valid_snps,
    const std::unordered_map<std::string, std::vector<uint64_t>>&
        gene_membership,
    const std::vector<double>& p_threshold, const size_t num_set)
{
    std::unordered_map<std::string, std::vector<uint64_t>> snps;
    bool is_gz = eqtl_name.substr(eqtl_name.find_last_of(".") + 1) == ("gz");
    std::ifstream eqtl;
    GZSTREAM_NAMESPACE::igzstream gz_eqtl;
    std::string error_message = "";
    std::string line;
    if (is_gz)
    {
        gz_eqtl.open(eqtl_name.c_str());
        if (!gz_eqtl.good())
        {
            error_message = "Error: Cannot open gziped eQTL file: " + eqtl_name;
            throw std::runtime_error(error_message);
        }
        std::getline(gz_eqtl, line);
    }
    else
    {
        eqtl.open(eqtl_name.c_str());
        if (!eqtl.is_open())
        {
            error_message = "Error: Cannot open eQTL file: " + eqtl_name;
            throw std::runtime_error(error_message);
        }
        std::getline(eqtl, line);
    }
    misc::trim(line);
    if (line.empty())
    {
        error_message = "Error: Header of eQTL file is empty";
        throw std::runtime_error(error_message);
    }

    size_t snp_idx = 0, p_idx = 0, gene_idx = 0;
    bool snp_found = false, p_found = false, gene_found = false;
    std::vector<std::string> token = misc::split(line);
    for (size_t i = 0; i < token.size(); ++i)
    {
        if (token[i] == eqtl_snp_id)
        {
            snp_idx = i;
            snp_found = true;
        }
        if (token[i] == eqtl_gene_id)
        {
            gene_idx = i;
            gene_found = true;
        }
        if (token[i] == eqtl_p)
        {
            p_idx = i;
            p_found = true;
        }
    }
    if (!snp_found || !p_found || !gene_found)
    {
        error_message =
            "Error: Required columns not found in the eQTL file: " + line;
        throw std::runtime_error(error_message);
    }
    size_t max_idx = 0;
    if (snp_idx > max_idx) max_idx = snp_idx;
    if (gene_idx > max_idx) max_idx = gene_idx;
    if (p_idx > max_idx) max_idx = p_idx;
    // now process the eqtl file
    size_t total_entry = 0, included_entry = 0;
    std::vector<std::string> snp_id;
    std::string cur_id, gene_name;
    double p_value = 2.0, pthres;
    // don't need to + 1 here, as background is included
    const size_t num_threshold = p_threshold.size();
    const size_t flag_size = ((num_set * num_threshold) / 64) + 1;
    std::vector<std::string> gene_name_info;
    while ((is_gz && getline(gz_eqtl, line)) || (!is_gz && getline(eqtl, line)))
    {
        misc::trim(line);
        if (line.empty()) continue;
        total_entry++;
        token = misc::split(line);
        if (max_idx >= token.size())
        {
            error_message = "Error: eQTL file does not contain enough columns!";
            throw std::runtime_error(error_message);
        }
        snp_id = misc::split(token[snp_idx], "_");
        // only use chr and bp to match
        if (snp_id.size() < 2)
        {
            error_message =
                "Error: Malform SNP ID. Should have format CHR_BP_A1_A2";
            throw std::runtime_error(error_message);
        }
        cur_id = snp_id[0] + "_" + snp_id[1];
        auto&& res = valid_snps.find(cur_id);
        if (res == valid_snps.end()) continue;
        // found
        included_entry++;
        gene_name = token[gene_idx];
        try
        {
            p_value = misc::convert<double>(token[p_idx]);
        }
        catch (const std::runtime_error&)
        {
            // check if it is NA
            std::transform(token[p_idx].begin(), token[p_idx].end(),
                           token[p_idx].begin(), ::toupper);
            if (token[p_idx] != "NAN" && token[p_idx] != "NA"
                && token[p_idx] != "NULL")
            {
                error_message =
                    "Error: Non-numeric p-value: " + misc::to_string(p_value);
                throw std::runtime_error(error_message);
            }
        }
        auto&& gene_info = gene_membership.find(gene_name);
        // Gene isn't found in MSigDB, so should be consider as part of the
        // background
        if (gene_info == gene_membership.end())
        {
            // This might be because of the new gene id format used by gtex
            // which follow the format geneID.Version
            gene_name_info = misc::split(gene_name, ".");
            gene_info = gene_membership.find(gene_name_info.front());
        }
        bool background_only = (gene_info == gene_membership.end());

        auto&& idx = snps.find(cur_id);
        if (idx == snps.end())
        {
            // num_set should contain the background here
            // use res->second for the RS ID instead of the eQTL ID
            snps[res->second] = std::vector<uint64_t>(flag_size, 0ULL);
        }
        // take care of background first
        size_t col = calculate_column(p_value, p_threshold, pthres);
        for (size_t col_idx = col; col_idx < num_threshold; ++col_idx)
        { snps[res->second][col_idx / 64] |= 1ULL << ((col_idx) % 64); }
        if (!background_only)
        {
            for (size_t set_idx = 1; set_idx < num_set; ++set_idx)
            {
                // loop through the remaining sets
                if ((gene_info->second[(set_idx) / 64] >> ((set_idx) % 64)) & 1)
                {
                    for (size_t col_idx = col; col_idx < num_threshold;
                         ++col_idx)
                    {
                        snps[res->second]
                            [(set_idx * num_threshold + col_idx) / 64] |=
                            1ULL << ((set_idx * num_threshold + col_idx) % 64);
                    }
                }
            }
        }
    }
    if (is_gz)
        gz_eqtl.close();
    else
        eqtl.close();
    // SNP matrix is way too big.
    // now we have SNPs, which contain the p-value of each SNP to each set
    std::cerr << snps.size() << " SNPs remaining after reading the eQTL file"
              << std::endl;
    return snps;
}

void generate_snp_sets(
    const std::unordered_map<std::string, std::vector<uint64_t>>& snps,
    const std::vector<std::string>& set_name,
    const std::vector<double>& p_thresholds, const std::string& out)
{
    std::ofstream output;
    output.open(out.c_str());
    std::string error_message = "";
    if (!output.is_open())
    {
        error_message = "Error: Cannot open output file - " + out + " to write";
        throw std::runtime_error(error_message);
    }
    // row = pathway, col = p-value threshold
    misc::vec2d<std::string> pathway_map(set_name.size(), p_thresholds.size(),
                                         "");
    std::cerr << "A total of " << set_name.size() * p_thresholds.size()
              << " SNP Sets should be generated" << std::endl;
    std::string rsid = "";
    const size_t num_sets = set_name.size();
    const size_t num_threshold = p_thresholds.size();
    for (auto& snp : snps)
    {
        // go through each SNP
        rsid = snp.first;
        for (size_t i_row = 0; i_row < num_sets; ++i_row)
        {
            for (size_t i_col = 0; i_col < num_threshold; ++i_col)
            {
                if ((snp.second[(i_row * num_threshold + i_col) / 64]
                     >> ((i_row * num_threshold + i_col) % 64))
                    & 1)
                { pathway_map(i_row, i_col).append(" " + rsid); }
            }
        }
    }
    // now we can write the output
    std::string cur_name;
    for (size_t row_idx = 0; row_idx < set_name.size(); ++row_idx)
    {
        cur_name = set_name[row_idx];
        std::replace(cur_name.begin(), cur_name.end(), '-', '_');
        for (size_t col_idx = 0; col_idx < p_thresholds.size(); ++col_idx)
        {
            // we don't need the space here, because pathway_map already
            // contain the space
            output << cur_name << "-" << p_thresholds[col_idx]
                   << pathway_map(row_idx, col_idx) << "\n";
        }
    }
    output.close();
}
#endif // FUNCTIONS_H
