#include "commander.h"

Commander::Commander() {}
Commander::~Commander() {}
bool Commander::initialize(int argc, char** argv)
{

    if (argc == 1)
    {
        usage();
        return false;
    }
    static const char* optString = "p:e:s:m:g:t:G:S:P:r:o:h?";
    static const struct option longOpts[] = {
        {"eqtl", required_argument, nullptr, 'e'},
        {"sumstat", required_argument, nullptr, 's'},
        {"msigdb", required_argument, nullptr, 'm'},
        {"gtf", required_argument, nullptr, 'g'},
        {"target", required_argument, nullptr, 't'},
        {"eSNP", required_argument, nullptr, 'S'},
        {"eGene", required_argument, nullptr, 'G'},
        {"eP", required_argument, nullptr, 'P'},
        {"rsid", required_argument, nullptr, 'r'},
        {"pthres", required_argument, nullptr, 'p'},
        {"out", required_argument, nullptr, 'o'},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}};

    int longIndex = 0;
    int opt = 0;
    opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    bool error = false;
    std::string error_message = "";
    while (opt != -1)
    {
        switch (opt)
        {
        case 'e': m_eqtl_name = optarg; break;
        case 's': m_sumstat_name = optarg; break;
        case 'm': m_msigdb_name = optarg; break;
        case 'g': m_gtf_name = optarg; break;
        case 't': m_target_name = optarg; break;
        case 'S': m_e_snp_id = optarg; break;
        case 'G': m_e_gene_id = optarg; break;
        case 'P': m_e_pvalue = optarg; break;
        case 'r': m_s_snp_id = optarg; break;
        case 'p':

        {
            // convert to vector double
            std::vector<std::string> token = misc::split(optarg);
            for (auto&& t : token)
            {
                try
                {
                    m_p_threshold.push_back(misc::convert<double>(t));
                }
                catch (...)
                {
                    error = true;
                    error_message.append(
                        "Error: Non-numeric p-value threshold observed\n");
                }
            }
        }
        break;
        case 'o': m_out = optarg; break;
        case 'h':
        case '?': usage(); exit(0);
        default:
            throw "Undefined operator, please use --help for more information!";
        }
        opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    }
    if (m_p_threshold.empty())
    { m_p_threshold = {5e-8, 5e-7, 5e-6, 5e-5, 5e-4, 0.005, 0.05, 0.5, 1}; }
    if (m_eqtl_name.empty())
    {
        error = true;
        error_message.append("Error: You must provide the eQTL file\n");
    }
    if (m_sumstat_name.empty())
    {
        error = true;
        error_message.append(
            "Error: You must provide the summary statistic file\n");
    }
    if (m_target_name.empty())
    {
        error = true;
        error_message.append("Error: You must provide the target file\n");
    }
    if (m_msigdb_name.empty())
    {
        error = true;
        error_message.append("Error: You must provide the MSigDB file\n");
    }
    if (m_gtf_name.empty())
    {
        error = true;
        error_message.append("Error: You must provide the GTF file\n");
    }
    if (m_e_snp_id.empty())
    {
        error = true;
        error_message.append(
            "Error: You must provide the SNP column name in the eQTL file\n");
    }
    if (m_e_gene_id.empty())
    {
        error = true;
        error_message.append(
            "Error: You must provide the Gene column name in the eQTL file\n");
    }
    if (m_e_pvalue.empty())
    {
        error = true;
        error_message.append("Error: You must provide the P-value column name "
                             "in the eQTL file\n");
    }
    if (m_s_snp_id.empty())
    {
        error = true;
        error_message.append("Error: You must provide the SNP column name in "
                             "the summary statistic file\n");
    }
    if (m_out.empty())
    {
        error = true;
        error_message.append("Error: You must provide the output file name\n");
    }
    if (error) { throw std::runtime_error(error_message); }
    std::cerr << "\neSNP-Sets " + m_version + " (" + m_date + ") \n";
    std::cerr << "https://github.com/choishingwan/eSNP-Sets\n";
    std::cerr << "(C) 2019 Shing Wan (Sam) Choi\n";
    std::cerr << "MIT License\n\n";
    std::cerr << argv[0] << " \\" << std::endl;
    std::cerr << "    --eGene " << m_e_gene_id << " \\" << std::endl;
    std::cerr << "    --eP " << m_e_pvalue << " \\" << std::endl;
    std::cerr << "    --eqtl " << m_eqtl_name << " \\" << std::endl;
    std::cerr << "    --eSNP " << m_e_snp_id << " \\" << std::endl;
    std::cerr << "    --gtf " << m_gtf_name << " \\" << std::endl;
    std::cerr << "    --msigdb " << m_msigdb_name << " \\" << std::endl;
    std::cerr << "    --out " << m_out << " \\" << std::endl;
    std::cerr << "    --pthres " << m_p_threshold.front();
    for (size_t i = 1; i < m_p_threshold.size(); ++i)
    { std::cerr << "," << m_p_threshold[i]; }
    std::cerr << std::endl;
    std::cerr << "    --rsid " << m_s_snp_id << " \\" << std::endl;
    std::cerr << "    --sumstat " << m_sumstat_name << " \\" << std::endl;
    std::cerr << "    --target " << m_target_name << std::endl;
    std::sort(m_p_threshold.begin(), m_p_threshold.end());
    return true;
}
void Commander::usage()
{
    std::cerr << "\neSNP-Sets " + m_version + " (" + m_date + ") \n";
    std::cerr << "https://github.com/choishingwan/eSNP-Sets\n";
    std::cerr << "(C) 2019 Shing Wan (Sam) Choi\n";
    std::cerr << "MIT License\n\n";
    std::cerr << "usage: eSNP-Sets [options] <-e eQTL file> <-t Target file> "
                 "<-s Sumstat file>"
              << std::endl;
    std::cerr << "    --eqtl     | -e    eQTL file" << std::endl;
    std::cerr << "    --sumstat  | -s    Summary statistic file" << std::endl;
    std::cerr << "    --msigdb   | -m    MSigDB file" << std::endl;
    std::cerr << "    --gtf      | -g    GTF file" << std::endl;
    std::cerr << "    --target   | -t    Target file" << std::endl;
    std::cerr << "    --eSNP     | -S    SNP ID column in the eQTL file"
              << std::endl;
    std::cerr << "    --eGene    | -G    Gene ID column in the eQTL file"
              << std::endl;
    std::cerr << "    --eP       | -P    P-value column in the eQTL file"
              << std::endl;
    std::cerr << "    --rsid     | -r    SNP ID column in the sumstat file"
              << std::endl;
    std::cerr << "    --pthres   | -p    p-value threshold to consider"
              << std::endl;
    std::cerr << "    --out      | -o    Output file name" << std::endl;
    std::cerr << "    --help     | -h    Display this message" << std::endl;
}
