#include <iostream>
#include <string>
#include <fstream>
#include <streambuf>

#include "denoising.h"
#include "point_cloud_io.h"

using namespace std;

void export_features(vector<FEATURE> & feature_classification, string filename)
{
	ofstream file;
	file.open(filename);
	for (int i = 0;i < feature_classification.size();i++)
	{
		if (feature_classification[i] == FLAT)
		{
			file << "0" << endl;
		}
		else
		{
			file << "1" << endl;
		}
	}
	file.close();
}

void read_parameters(DenoisingParameters & parameters, string csv_filename)
{
	std::ifstream fin(csv_filename);
	vector<string> row;
	string line, word, temp;

	while (getline(fin, line))
	{
		row.clear();

		stringstream s(line);
		while (getline(s, word, ',')) 
		{
			row.push_back(word);
		}
		if (row.size() == 2)
		{
			string param = row[0];
			string value = row[1];
			
			if (param == "n_ext")
				parameters.n_ext = stoi(value);
			else if (param == "n_int")
				parameters.n_int = stoi(value);
			else if (param == "k")
				parameters.k = stoi(value);
			else if (param == "r_r")
				parameters.r_r = stod(value);
			else if (param == "alpha")
				parameters.alpha = stod(value);
			else if (param == "beta")
				parameters.beta = stod(value);
			else if (param == "gamma")
				parameters.gamma = stod(value);
			else if (param == "a_0_ratio")
				parameters.a_0_ratio = stod(value);
			else if (param == "tau_n")
				parameters.tau_n = stod(value);
			else if (param == "n_ns")
				parameters.n_ns = stoi(value);
			else if (param == "sigma_nn")
				parameters.sigma_nn = stod(value);
			else if (param == "sigma_ns_ratio")
				parameters.sigma_ns_ratio = stod(value);
			else if (param == "n_nc")
				parameters.n_nc = stoi(value);
			else if (param == "delta_cc_ratio")
				parameters.delta_cc_ratio = stod(value);
			else if (param == "theta_half")
				parameters.theta_half = stod(value);
			else if (param == "n_fp")
				parameters.n_fp = stoi(value);
			else if (param == "upsilon_f")
				parameters.upsilon_f = stod(value);
			else if (param == "upsilon_ec")
				parameters.upsilon_ec = stod(value);
			else if (param == "tau_o_ratio")
				parameters.tau_o_ratio = stod(value);
			else
				cout << "Unexpected parameter: " << param << endl;
		}		
	}
}


int main(int argc, char **argv)
{
	if (argc == 3)
	{
		std::string input_file_name(argv[1]);
		std::string output_file_name(argv[2]);
		cout << input_file_name << endl << output_file_name << endl;
		PointCloud noisy_pc;
		import_point_cloud(noisy_pc, input_file_name);
		DenoisingParameters denoising_parameters;
		DenoisingData denoising_data;
		denoising_parameters.show();
		PointCloud denoised_pc = denoise(noisy_pc, denoising_parameters, denoising_data);
		export_point_cloud(denoised_pc, output_file_name);
		return 0;
	}
	else if (argc == 4)
	{
		std::string input_file_name(argv[1]);
		std::string output_file_name(argv[2]);
		std::string params_file_name(argv[3]);
		cout << input_file_name << endl << output_file_name << endl;
		PointCloud noisy_pc;
		import_point_cloud(noisy_pc, input_file_name);
		DenoisingParameters denoising_parameters;
		read_parameters(denoising_parameters, params_file_name);
		DenoisingData denoising_data;
		denoising_parameters.show();
		PointCloud denoised_pc = denoise(noisy_pc, denoising_parameters, denoising_data);
		export_point_cloud(denoised_pc, output_file_name);
		return 0;
	}
	else if (argc == 5)
	{
		std::string input_file_name(argv[1]);
		std::string output_file_name(argv[2]);
		std::string params_file_name(argv[3]);
		std::string features_file_name(argv[4]);
		cout << input_file_name << endl << output_file_name << endl;
		PointCloud noisy_pc;
		import_point_cloud(noisy_pc, input_file_name);
		DenoisingParameters denoising_parameters;
		read_parameters(denoising_parameters, params_file_name);
		DenoisingData denoising_data;
		denoising_parameters.show();
		PointCloud denoised_pc = denoise(noisy_pc, denoising_parameters, denoising_data);
		export_point_cloud(denoised_pc, output_file_name);
		export_features(denoising_data.feature_classification, features_file_name);
		return 0;
	}
	else
	{
		std::cout << "Usage:" << std::endl;
		std::cout << "\t pcdenoising point_cloud.ply point_cloud_denoised.ply" << std::endl;
		std::cout << "\t pcdenoising point_cloud.ply point_cloud_denoised.ply point_cloud_parameters.csv" << std::endl;
		std::cout << "\t pcdenoising point_cloud.ply point_cloud_denoised.ply point_cloud_parameters.csv point_cloud_features.txt" << std::endl;
		return 1;
	}
}