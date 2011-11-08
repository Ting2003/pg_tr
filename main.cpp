#include "main.h"

const char * usage="Usage: %s [-eorbILifl] benchmark\n\
    -e EPSILON\n\
    -o OMEGA\n\
    -r overlap ratio\n\
    -b max block nodes\n\
    -I block iterative (default)\n\
    -L direct LU\n\
    -i input file\n\
    -f output file\n\
    -l log file (default to screen)\n"
;

const char * usage2="Usage: %s -i input -f output\n";

int main(int argc, char * argv[]){
#if 0
	int c;
	int mode=0;
	double epsilon, omega, overlap_ratio;
	size_t max_block_nodes;
	//char * logfile="/dev/null";
	char * logfile=NULL;
	char * input=NULL, * output=NULL;
	bool input_flag = false, output_flag = false;
	Circuit::get_parameters(epsilon, omega, overlap_ratio, 
			max_block_nodes, mode);

	while( ( c = getopt(argc, argv, "i:f:e:o:r:b:l:LI")) != -1 ){
		switch(c){
		case 'e':
			epsilon = atof(optarg);
			break;
		case 'o':
			omega = atof(optarg);
			break;
		case 'r':
			overlap_ratio = atof(optarg);
			break;
		case 'b':
			max_block_nodes = atof(optarg);
			break;
		case 'L':
			mode = 1;
			break;
		case 'I':
			mode = 0;
			break;
		case 'l':
			logfile = optarg;
			break;
		case 'i':
			input = optarg;
			input_flag = true;
			break;
		case 'f':
			output = optarg;
			output_flag = true;
			break;
		case '?':
		default:
			fprintf(stderr, usage2, argv[0]);
			exit(EXIT_FAILURE);
		}
	}
	//if( argc == optind ) report_exit(usage2);
	if( !input_flag || ! output_flag ){
		fprintf(stderr, usage2, argv[0]);
		exit(EXIT_FAILURE);
	}

	open_logfile(logfile);
	if( freopen(output, "w", stdout) == NULL )
		report_exit("Ouptut file error\n");

	Circuit::set_parameters(epsilon, omega, overlap_ratio, 
			max_block_nodes, mode);

	Tran tran;
	clock_t s, e;
	s = clock();
	// start to parfile
	vector<Circuit *> cktlist;
	Parser parser(&cktlist);
	parser.parse(input, tran);
	e = clock();
	clog<<"parse cost, omp: "<<1.0*(e-s)/CLOCKS_PER_SEC<<endl;
	return 0;
	// do the job
	for(size_t i=0;i<cktlist.size();i++){
		Circuit * ckt = cktlist[i];
		ckt->solve(tran);
	        // after that, this circuit can be released
		delete ckt;
	}
	//tran.print_tr_nodes();

	close_logfile();

#endif
	clock_t t1, t2;
	t1 = clock();
	size_t a = 0;
	for(size_t i=0;i<1e10;i++)
		a = i;
	t2 = clock();
	cout<<"simple cost: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
	return 0;
}
