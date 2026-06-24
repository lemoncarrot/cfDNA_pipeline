//fast calculation of p-values in c++
//author: billy lau
//
//compile me with g++ -Ofast -std=c++1z -o cpg_level_ttest cpg_level_ttest.cpp
//
//Usage: CREST_pvalues_only_cpp [input methylation file] [masking definitions] [output file] [min samples per group]
//  min samples per group: a CpG is only tested if it has >= this many non-NA samples in
//  BOTH the case and the control group; CpGs below the threshold are skipped (not output).
//test with
//nohup ./CREST_pvalues_only_cpp /mnt/ix1/Projects/M081_20210811_CREST_proc/A01_CREST_sequencing_data/bedgraph_filter/CREST.aggr.union.filter.bedgraph test_samples_oob.txt test_pvalues.txt 5

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <boost/math/distributions/students_t.hpp>
#include <vector>
#include <numeric>
#include <chrono>
#include <ctime>
#include <unordered_map>

double calculateVariance(std::vector<double> &data, double mean) {
    double sum = 0;
    for (double x : data) {
        sum += (x - mean) * (x - mean);
    }
    return sum / (data.size() - 1);
}

double calculateTStatistic(double mean1, double mean2, double var1, double var2, int size1, int size2) {
    return (mean1 - mean2)/sqrt(var1/size1 + var2/size2);
}

double calculatePValue(double tStatistic, int size1, int size2, double var1, double var2) {
    using boost::math::cdf;
    using boost::math::complement;
    using boost::math::students_t;
    
    double df = (var1/size1 + var2/size2) * (var1/size1 + var2/size2) / ((var1 * var1)/(size1 * size1 * (size1 - 1)) + (var2 * var2)/(size2 * size2 * (size2 - 1)));
    students_t dist(df);
    return 2 * cdf(complement(dist, abs(tStatistic))); //two-sided t-test
}

int main(int argc, char* argv[]) {
    std::cin.sync_with_stdio(false);

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " [input methylation file] [masking definitions] [output file]\n";
        return 1;
    }

    std::string line;
    std::ifstream file_meth;
    file_meth.open(argv[1]);
    std::ifstream file_groups;
    file_groups.open(argv[2]);
    std::ofstream file_pvalues;
    file_pvalues.open(argv[3]);
    
    // build name -> label map from sample_def.txt
    std::unordered_map<std::string, std::string> group_map;
    while (getline(file_groups, line)) {
    	std::stringstream ss(line);
    	std::string name, label;
    	std::getline(ss, name, '\t');
    	std::getline(ss, label, '\t');
    	group_map[name] = label;
    }
    file_groups.close();

    // parse header row into col_names for name-based lookup
    std::vector<std::string> col_names;
    std::getline(file_meth, line);
    {
        std::stringstream ss(line);
        std::string col;
        while (std::getline(ss, col, '\t')) col_names.push_back(col);
    }

    std::vector<double> case_group, control_group;
    int line_num = 1;

    //https://en.cppreference.com/w/cpp/chrono
    auto cur_time = std::chrono::system_clock::now();

    file_pvalues << "chrom\tstart\tend\tt-statistic\tpval\tmean_case\tmean_control\tn_case\tn_control" << "\n";

    while (std::getline(file_meth, line)) {
        std::stringstream ss(line);
        std::string item;

        case_group.clear();
        control_group.clear();

        //output benchmarking time
        if (line_num % 100000 == 0) {
        	auto end_time = std::chrono::system_clock::now();
        	std::chrono::duration<double> elapsed_seconds = end_time - cur_time;
        	std::cerr << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";
        	cur_time = std::chrono::system_clock::now();
        }

        //populate the cases and controls
        int index = 0;
        std::string in_chr,in_start,in_end;
        while (std::getline(ss, item, '\t')) {
        	if (index == 0) {in_chr = item;} else
        	if (index == 1) {in_start = item;} else
        	if (index == 2) {in_end = item;} else
        	if (group_map[col_names[index]] == "mask") {} else
        	if (group_map[col_names[index]] == "oob") {} else
        	if (item == "NA") {} else
        	if (group_map[col_names[index]] == "case") {case_group.push_back(std::stod(item));} else
        	if (group_map[col_names[index]] == "control") {control_group.push_back(std::stod(item));}
            index++;
        }

        double mean_case = std::reduce(case_group.begin(), case_group.end()) / case_group.size();
        double mean_control = std::reduce(control_group.begin(), control_group.end()) / control_group.size();

        double var_case = calculateVariance(case_group, mean_case);
        double var_control = calculateVariance(control_group, mean_control);
        int size_case = case_group.size();
        int size_control = control_group.size();
        
        //take care of edge cases to speed up
        line_num += 1;
        if ((size_case <= 1) || (size_control <= 1)) { //not enough samples
        	file_pvalues << in_chr << "\t" << in_start << "\t" << in_end << "\t" << "--" << "\t" << "--" << "\t" << mean_case << "\t" << mean_control << "\t" << size_case << "\t" << size_control << "\n";
        } else if ((mean_case - mean_control == 0) && (var_case == 0) && (var_control == 0)) { //no diff in variance and mean; nan p value
        	file_pvalues << in_chr << "\t" << in_start << "\t" << in_end << "\t" << "nan" << "\t" << "nan" << "\t" << mean_case << "\t" << mean_control << "\t" << size_case << "\t" << size_control << "\n";
        } else if ((var_case == 0) && (var_control == 0)) { //no diff in variance, different mean; -inf t-statistic and 0 p-value
        	file_pvalues << in_chr << "\t" << in_start << "\t" << in_end << "\t" << "-inf" << "\t" << "0" << "\t" << mean_case << "\t" << mean_control << "\t" << size_case << "\t" << size_control << "\n";
        } else { //finally calculate p value
			double tStatistic = calculateTStatistic(mean_case, mean_control, var_case, var_control, size_case, size_control);
			double pValue = calculatePValue(tStatistic, size_case, size_control, var_case, var_control);
		
			file_pvalues << in_chr << "\t" << in_start << "\t" << in_end << "\t" << tStatistic << "\t" << pValue << "\t" << mean_case << "\t" << mean_control << "\t" << size_case << "\t" << size_control << "\n";
        }
        
    }
    file_meth.close();
    file_pvalues.close();
    return 0;
}
