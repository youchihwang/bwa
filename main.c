#include <stdio.h>
#include <string.h>
#include "kstring.h"
#include "utils.h"

// <<< youchihwang, 2019/01/16, computing left right length from seed before smith waterman ksw_extend2().
#include <inttypes.h>
#include <time.h>
// >>> youchihwang, 2019/01/16, computing left right length from seed before smith waterman ksw_extend2().

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.7.17-r1194-dirty"
#endif

int bwa_fa2pac(int argc, char *argv[]);
int bwa_pac2bwt(int argc, char *argv[]);
int bwa_bwtupdate(int argc, char *argv[]);
int bwa_bwt2sa(int argc, char *argv[]);
int bwa_index(int argc, char *argv[]);
int bwt_bwtgen_main(int argc, char *argv[]);

int bwa_aln(int argc, char *argv[]);
int bwa_sai2sam_se(int argc, char *argv[]);
int bwa_sai2sam_pe(int argc, char *argv[]);

int bwa_bwtsw2(int argc, char *argv[]);

int main_fastmap(int argc, char *argv[]);
int main_mem(int argc, char *argv[]);
int main_shm(int argc, char *argv[]);

int main_pemerge(int argc, char *argv[]);
int main_maxk(int argc, char *argv[]);
	
static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: bwa (alignment via Burrows-Wheeler transformation)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Heng Li <lh3@sanger.ac.uk>\n\n");
	fprintf(stderr, "Usage:   bwa <command> [options]\n\n");
	fprintf(stderr, "Command: index         index sequences in the FASTA format\n");
	fprintf(stderr, "         mem           BWA-MEM algorithm\n");
	fprintf(stderr, "         fastmap       identify super-maximal exact matches\n");
	fprintf(stderr, "         pemerge       merge overlapping paired ends (EXPERIMENTAL)\n");
	fprintf(stderr, "         aln           gapped/ungapped alignment\n");
	fprintf(stderr, "         samse         generate alignment (single ended)\n");
	fprintf(stderr, "         sampe         generate alignment (paired ended)\n");
	fprintf(stderr, "         bwasw         BWA-SW for long queries\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "         shm           manage indices in shared memory\n");
	fprintf(stderr, "         fa2pac        convert FASTA to PAC format\n");
	fprintf(stderr, "         pac2bwt       generate BWT from PAC\n");
	fprintf(stderr, "         pac2bwtgen    alternative algorithm for generating BWT\n");
	fprintf(stderr, "         bwtupdate     update .bwt to the new format\n");
	fprintf(stderr, "         bwt2sa        generate SA from BWT and Occ\n");
	fprintf(stderr, "\n");
	fprintf(stderr,
"Note: To use BWA, you need to first index the genome with `bwa index'.\n"
"      There are three alignment algorithms in BWA: `mem', `bwasw', and\n"
"      `aln/samse/sampe'. If you are not sure which to use, try `bwa mem'\n"
"      first. Please `man ./bwa.1' for the manual.\n\n");
	return 1;
}

// <<< youchihwang, 2019/01/16, computing left right length from seed before smith waterman ksw_extend2().
extern int extern_n_seqs;
extern int64_t total_left_length;
extern int64_t total_right_length;
extern int64_t max_left_length;
extern int64_t min_left_length;
extern int64_t max_right_length;
extern int64_t min_right_length;
// >>> youchihwang, 2019/01/16, computing left right length from seed before smith waterman ksw_extend2().

int main(int argc, char *argv[])
{
	extern char *bwa_pg;
	int i, ret;
	double t_real;
	kstring_t pg = {0,0,0};
	t_real = realtime();
	ksprintf(&pg, "@PG\tID:bwa\tPN:bwa\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
	for (i = 1; i < argc; ++i) ksprintf(&pg, " %s", argv[i]);
	bwa_pg = pg.s;
	if (argc < 2) return usage();
	if (strcmp(argv[1], "fa2pac") == 0) ret = bwa_fa2pac(argc-1, argv+1);
	else if (strcmp(argv[1], "pac2bwt") == 0) ret = bwa_pac2bwt(argc-1, argv+1);
	else if (strcmp(argv[1], "pac2bwtgen") == 0) ret = bwt_bwtgen_main(argc-1, argv+1);
	else if (strcmp(argv[1], "bwtupdate") == 0) ret = bwa_bwtupdate(argc-1, argv+1);
	else if (strcmp(argv[1], "bwt2sa") == 0) ret = bwa_bwt2sa(argc-1, argv+1);
	else if (strcmp(argv[1], "index") == 0) ret = bwa_index(argc-1, argv+1);
	else if (strcmp(argv[1], "aln") == 0) ret = bwa_aln(argc-1, argv+1);
	else if (strcmp(argv[1], "samse") == 0) ret = bwa_sai2sam_se(argc-1, argv+1);
	else if (strcmp(argv[1], "sampe") == 0) ret = bwa_sai2sam_pe(argc-1, argv+1);
	else if (strcmp(argv[1], "bwtsw2") == 0) ret = bwa_bwtsw2(argc-1, argv+1);
	else if (strcmp(argv[1], "dbwtsw") == 0) ret = bwa_bwtsw2(argc-1, argv+1);
	else if (strcmp(argv[1], "bwasw") == 0) ret = bwa_bwtsw2(argc-1, argv+1);
	else if (strcmp(argv[1], "fastmap") == 0) ret = main_fastmap(argc-1, argv+1);
	else if (strcmp(argv[1], "mem") == 0) ret = main_mem(argc-1, argv+1);
	else if (strcmp(argv[1], "shm") == 0) ret = main_shm(argc-1, argv+1);
	else if (strcmp(argv[1], "pemerge") == 0) ret = main_pemerge(argc-1, argv+1);
	else if (strcmp(argv[1], "maxk") == 0) ret = main_maxk(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}

    // <<< youchihwang, 2019/01/16, computing left right length from seed before smith waterman ksw_extend2().
    FILE * fp_w = fopen("bohao", "a+w");
    time_t nowtime;
    nowtime = time(NULL);
    char tmp[64];
    strftime(tmp, sizeof(tmp), "%Y-%m-%d %H %M %S", localtime(&nowtime));

    if (fp_w == NULL)
        return -1;

    // printf("\n");
    // printf("\n");
    // printf("min left length : %" PRId64 "\n", min_left_length);
    // printf("max left length : %" PRId64 "\n", max_left_length);
    // printf("min right length : %" PRId64 "\n", min_right_length);
    // printf("max right length : %" PRId64 "\n", max_right_length);
    // printf("The total read : %d \n", extern_n_seqs);
    // printf("The total left length : %" PRId64 "\n", total_left_length);
    // printf("The total right length : %" PRId64 "\n", total_right_length);
    // printf("The average left length of each read : %" PRId64 "\n", total_left_length/extern_n_seqs);
    // printf("The average right length of each read : %" PRId64 "\n", total_right_length/extern_n_seqs);
    // printf("\n");
    // printf("\n");

    fprintf(fp_w, "\n");
    fprintf(fp_w, "date time : %s \n", tmp);
    fprintf(fp_w, "min left length : %" PRId64 "\n", min_left_length);
    fprintf(fp_w, "max left length : %" PRId64 "\n", max_left_length);
    fprintf(fp_w, "min right length : %" PRId64 "\n", min_right_length);
    fprintf(fp_w, "max right length : %" PRId64 "\n", max_right_length);
    fprintf(fp_w, "The total read : %d \n", extern_n_seqs);
    fprintf(fp_w, "The total left length : %" PRId64 "\n", total_left_length);
    fprintf(fp_w, "The total right length : %" PRId64 "\n", total_right_length);
    fprintf(fp_w, "The average left length of each read : %" PRId64 "\n", total_left_length/extern_n_seqs);
    fprintf(fp_w, "The average right length of each read : %" PRId64 "\n", total_right_length/extern_n_seqs);
    fprintf(fp_w, "\n");

    fclose(fp_w);
    // >>> youchihwang, 2019/01/16, computing left right length from seed before smith waterman ksw_extend2().

	err_fflush(stdout);
	err_fclose(stdout);
	if (ret == 0) {
		fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
		fprintf(stderr, "[%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());
	}
	free(bwa_pg);
	return ret;
}
