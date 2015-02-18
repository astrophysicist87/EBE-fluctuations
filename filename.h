string outputfilename(int order_harm);
string outputfilename();

string outputfilename(int order_harm)
{
	stringstream out;
	out << order_harm;
	string order_string = out.str();
	out.str(std::string());

	time_t now = time(0);;
	tm *ltm = localtime(&now);

	out << 1900 + ltm->tm_year;
	string year_string = out.str();
	out.str(std::string());

	out << 1 + ltm->tm_mon;
	string month_string = out.str();
	out.str(std::string());

	out << ltm->tm_mday;
	string day_string = out.str();
	out.str(std::string());

	out << 1 + ltm->tm_hour;
	string hour_string = out.str();
	out.str(std::string());

	out << 1 + ltm->tm_min;
	string minute_string = out.str();
	out.str(std::string());

	out << 1 + ltm->tm_sec;
	string second_string = out.str();
	out.str(std::string());

	string output_filename = "R2_s_order_" + order_string + "_" + year_string + "_" + month_string + "_"
					+ day_string + "_" + hour_string + "_" + minute_string + "_" + second_string +"_4d.dat";

	return (output_filename);
}

string outputfilename()
{
	stringstream out;

	time_t now = time(0);;
	tm *ltm = localtime(&now);

	out << 1900 + ltm->tm_year;
	string year_string = out.str();
	out.str(std::string());

	out << 1 + ltm->tm_mon;
	string month_string = out.str();
	out.str(std::string());

	out << ltm->tm_mday;
	string day_string = out.str();
	out.str(std::string());

	out << 1 + ltm->tm_hour;
	string hour_string = out.str();
	out.str(std::string());

	out << 1 + ltm->tm_min;
	string minute_string = out.str();
	out.str(std::string());

	out << 1 + ltm->tm_sec;
	string second_string = out.str();
	out.str(std::string());

	string output_filename = "R2_s_order_" + year_string + "_" + month_string + "_"
					+ day_string + "_" + hour_string + "_" + minute_string + "_" + second_string +"_4d.dat";

	return (output_filename);
}
