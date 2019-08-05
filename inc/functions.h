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
void set_bit(const size_t& idx, std::vector<size_t>& flag)
{
    if (flag.empty())
    { throw std::runtime_error("Error: Empty Flag detected"); }
    else if (!(flag.size() > idx / BITCT))
    {
        std::string error_message =
            "Error: Flag index out of range: " + std::to_string(idx)
            + ". Max idx allowed is " + std::to_string(flag.size() * BITCT - 1);
        throw std::out_of_range(error_message);
    }
    flag[idx / BITCT] |= ONEUL << ((idx) % BITCT);
}

void range_set_bit(const size_t& idx, const size_t& end,
                   std::vector<size_t>& flag)
{
    if (flag.empty())
    { throw std::runtime_error("Error: Empty Flag detected"); }
    else if (!(flag.size() > idx / BITCT))
    {
        std::string error_message =
            "Error: Flag index out of range: " + std::to_string(idx)
            + ". Max idx allowed is " + std::to_string(flag.size() * BITCT - 1);
        throw std::out_of_range(error_message);
    }
    else if (flag.size() * BITCT < end)
    {
        std::string error_message =
            "Error: Flag end bound out of range: " + std::to_string(end)
            + ". Max idx allowed is " + std::to_string(flag.size() * BITCT);
        throw std::out_of_range(error_message);
    }
    for (size_t start = idx; start < end; ++start) { set_bit(start, flag); }
}

bool get_bit(const size_t& idx, const std::vector<size_t>& flag)
{
    assert(flag.size() > idx / BITCT);
    return ((flag[(idx) / BITCT] >> ((idx) % BITCT)) & 1);
}

bool is_na(const std::string& in)
{
    std::string input = in;
    std::transform(input.begin(), input.end(), input.begin(), ::toupper);
    return (input == "NAN" || input == "NA" || input == "NULL");
}

bool parse_attribute(const std::string& attribute_str, std::string& gene_id,
                     std::string& gene_name)
{
    assert(!attribute_str.empty());
    std::vector<std::string> attributes = misc::split(attribute_str, ";");
    std::vector<std::string> token;
    bool found_id = false, found_name = false;
    for (auto&& a : attributes)
    {
        // remove space before the attribute name
        misc::trim(a);
        if (a.rfind("gene_id", 0) == 0)
        {
            token = misc::split(a, " ");
            if (token.size() != 2)
            {
                throw std::runtime_error("Error: Malformed attribute value: "
                                         + a);
            }
            gene_id = token.back();
            gene_id.erase(std::remove(gene_id.begin(), gene_id.end(), '\"'),
                          gene_id.end());
            if (found_name) return true;
            found_id = true;
        }
        else if (a.rfind("gene_name", 0) == 0)
        {
            token = misc::split(a, " ");
            if (token.size() != 2)
            {
                throw std::runtime_error("Error: Malformed attribute value: "
                                         + a);
            }
            gene_name = token.back();
            gene_name.erase(
                std::remove(gene_name.begin(), gene_name.end(), '\"'),
                gene_name.end());
            if (found_id) return true;
            found_name = true;
        }
    }
    return false;
}

double get_p(const std::string& p)
{
    double p_value = 2.0;
    try
    {
        p_value = misc::convert<double>(p);
    }
    catch (const std::runtime_error&)
    {
        // check if it is NA
        if (!is_na(p))
        {
            std::string error_message = "Error: Non-numeric p-value: " + p;
            throw std::runtime_error(error_message);
        }
    }
    return p_value;
}

bool is_gz_file(const std::string& name)
{
    const unsigned char gz_magic[2] = {0x1f, 0x8b};
    FILE* fp;
    if ((fp = fopen(name.c_str(), "rb")) == nullptr)
    { throw std::runtime_error("Error: Cannot open file - " + name); }
    unsigned char buf[2];
    if (fread(buf, 1, 2, fp) == 2)
    {
        if (buf[0] == gz_magic[0] && buf[1] == gz_magic[1]) { return true; }
        return false;
    }
    else
    {
        // can open the file, but can't read the magic number.
        return false;
    }
}

bool open_file(const std::string& name, GZSTREAM_NAMESPACE::igzstream& gz,
               std::ifstream& file)
{
    bool is_gz = is_gz_file(name);
    if (is_gz)
    {
        gz.open(name.c_str());
        if (!gz.good())
        {
            std::string error_message =
                "Error: Cannot open gziped file - " + name;
            throw std::runtime_error(error_message);
        }
    }
    else
    {
        file.open(name.c_str());
        if (!file.is_open())
        {
            std::string error_message = "Error: Cannot open file - " + name;
            throw std::runtime_error(error_message);
        }
    }
    return is_gz;
}

bool get_idx(const std::vector<std::string>& token, const std::string& id,
             size_t& idx)
{
    idx = 0;
    for (; idx < token.size(); ++idx)
    {
        if (token[idx] == id) { return true; }
    }
    return false;
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

std::string get_bim_name(const std::string& bim)
{
    if (bim.find(".") != std::string::npos
        && bim.substr(bim.find_last_of(".") + 1) == ("bim"))
    { return bim; }
    else
    {
        return bim + ".bim";
    }
}
std::unordered_map<std::string, std::unordered_set<std::string>>
gen_gene_id_map(const std::string& gtf_name)
{
    std::unordered_map<std::string, std::unordered_set<std::string>> result;

    std::ifstream gtf;
    GZSTREAM_NAMESPACE::igzstream gz_gtf;
    std::string error_message;
    const bool is_gz = open_file(gtf_name, gz_gtf, gtf);
    std::string line;
    std::vector<std::string> token;
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
        if (parse_attribute(token[8], gene_id, gene_name))
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

std::unordered_map<std::string, std::vector<size_t>> gen_gene_membership(
    const std::string& msigdb_name,
    const std::unordered_map<std::string, std::unordered_set<std::string>>&
        gene_map,
    std::vector<std::string>& set_name)
{
    std::unordered_map<std::string, std::vector<size_t>> result;
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
                { set_bit(idx, result[gene_id]); }
                else
                {
                    std::vector<size_t> tmp(flag_size, 0ULL);
                    result[gene_id] = tmp;
                    set_bit(0, result[gene_id]);
                    set_bit(idx, result[gene_id]);
                    result[gene_id] = tmp;
                }
            }
        }
        ++idx;
    }
    msigdb.close();
    std::cerr << num_sets << " Set(s) found in the MSigDB file" << std::endl;
    std::cerr << "With total of " << result.size() << " genes" << std::endl;
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
    // bim file format
    const bool is_gz = open_file(sumstat_name, gz_sumstat, sumstat);
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
    if (!get_idx(token, snp_id, snp_idx))
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
    std::string bim_name = get_bim_name(target);
    bim.open(bim_name.c_str());
    if (!bim.is_open())
    {
        error_message = "Error: Cannot open bim file: " + bim_name;
        throw std::runtime_error(error_message);
    }
    size_t num_not_found = 0;
    size_t num_line = 0;
    const size_t chr_idx = 0, rs_idx = 1, bp_idx = 2;
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
        if (snp_in_sumstat.find(token[rs_idx]) != snp_in_sumstat.end())
        {
            // this snp is found in the summary statistic
            // we assume the chromosome ID doesn't starts with CHR
            std::string id = token[chr_idx] + "_" + token[bp_idx];
            if (result.find(id) != result.end())
            {
                std::cerr << "Warning: Duplicated SNP ID generated - " << id
                          << std::endl;
            }
            // now store the rsid (token[1]) with the eqtl id
            result[id] = token[rs_idx];
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


std::unordered_map<std::string, std::vector<size_t>> gen_binary_pathway_member(
    const std::string& eqtl_name, const std::string& eqtl_snp_id,
    const std::string& eqtl_gene_id, const std::string& eqtl_p,
    const std::unordered_map<std::string, std::string>& valid_snps,
    const std::unordered_map<std::string, std::vector<size_t>>& gene_membership,
    const std::vector<double>& p_threshold, const size_t num_set)
{
    std::unordered_map<std::string, std::vector<size_t>> snps;
    std::unordered_map<std::string, std::string>::const_iterator id_map;
    std::ifstream eqtl;
    GZSTREAM_NAMESPACE::igzstream gz_eqtl;
    const bool is_gz = open_file(eqtl_name, gz_eqtl, eqtl);
    std::string error_message = "";
    std::string line;
    if (is_gz) { std::getline(gz_eqtl, line); }
    else
    {
        std::getline(eqtl, line);
    }
    misc::trim(line);
    if (line.empty())
    {
        error_message = "Error: Header of eQTL file is empty";
        throw std::runtime_error(error_message);
    }
    size_t snp_idx = 0, p_idx = 0, gene_idx = 0;
    std::vector<std::string> token = misc::split(line);
    bool snp_found = get_idx(token, eqtl_snp_id, snp_idx);
    bool p_found = get_idx(token, eqtl_p, p_idx);
    bool gene_found = get_idx(token, eqtl_gene_id, gene_idx);
    if (!snp_found || !p_found || !gene_found)
    {
        error_message =
            "Error: Required columns not found in the eQTL file: " + line;
        throw std::runtime_error(error_message);
    }
    size_t max_idx = snp_idx;
    if (gene_idx > max_idx) max_idx = gene_idx;
    if (p_idx > max_idx) max_idx = p_idx;
    // now process the eqtl file
    // don't need to + 1 here, as background is included
    const size_t num_threshold = p_threshold.size();
    const size_t flag_size = ((num_set * num_threshold) / 64) + 1;
    size_t total_entry = 0, included_entry = 0;
    std::string eqtl_id, gene_name, rsid;
    double p_value = 2.0, pthres = 0.0;
    std::vector<std::string> snp_id;
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
        eqtl_id = snp_id[0] + "_" + snp_id[1];
        id_map = valid_snps.find(eqtl_id);
        if (id_map == valid_snps.end()) continue;
        rsid = id_map->second;
        // found
        included_entry++;
        gene_name = token[gene_idx];
        p_value = get_p(token[p_idx]);
        auto&& gene_info = gene_membership.find(gene_name);
        if (gene_info == gene_membership.end())
        {
            // This might be because of the new gene id format used by gtex
            // which follow the format geneID.Version
            gene_name_info = misc::split(gene_name, ".");
            gene_info = gene_membership.find(gene_name_info.front());
        }
        // Because we are handling eQTL data set, any SNP loaded from the eQTL
        // are "linked" to a gene, and therefore should be used as background
        bool background_only = (gene_info == gene_membership.end());

        auto&& snp_flag = snps.find(rsid);
        if (snp_flag == snps.end())
        {
            // num_set should contain the background here
            snps[rsid] = std::vector<size_t>(flag_size, 0ULL);
            snp_flag = snps.find(rsid);
        }
        // take care of background first
        size_t col = calculate_column(p_value, p_threshold, pthres);
        range_set_bit(col, num_threshold, snp_flag->second);
        if (!background_only)
        {
            // start at 1 because background has already been handled
            for (size_t set_idx = 1; set_idx < num_set; ++set_idx)
            {
                // loop through the remaining sets
                if (!get_bit(set_idx, gene_info->second)) continue;
                range_set_bit(set_idx * num_threshold + col,
                              set_idx * num_threshold + num_threshold,
                              snp_flag->second);
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
    const std::unordered_map<std::string, std::vector<size_t>>& snps,
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
                if (get_bit(i_row * num_threshold + i_col, snp.second))
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
        std::replace(cur_name.begin(), cur_name.end(), ' ', '_');
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
