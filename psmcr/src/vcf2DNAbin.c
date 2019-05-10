#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>

/* translation table CHAR -> DNAbin (from ape) */
static unsigned char tab_trans[] = {
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 0-9 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 10-19 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 20-29 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 30-39 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x04, 0x00, 0x00, 0x00, 0x00, /* 40-49 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 50-59 */
	0x00, 0x00, 0x00, 0x02, 0x00, 0x88, 0x70, 0x28, 0xd0, 0x00, /* 60-69 */
	0x00, 0x48, 0xb0, 0x00, 0x00, 0x50, 0x00, 0xa0, 0xf0, 0x00, /* 70-79 */
	0x00, 0x00, 0xc0, 0x60, 0x18, 0x00, 0xe0, 0x90, 0x00, 0x30, /* 80-89 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x88, 0x70, 0x28, /* 90-99 */
	0xd0, 0x00, 0x00, 0x48, 0xb0, 0x00, 0x00, 0x50, 0x00, 0xa0, /* 100-109 */
	0xf0, 0x00, 0x00, 0x00, 0xc0, 0x60, 0x18, 0x00, 0xe0, 0x90, /* 110-119 */
	0x00, 0x30, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 120-129 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 130-139 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 140-149 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 150-159 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 160-169 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 170-179 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 180-189 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 190-199 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 200-209 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 210-219 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 220-229 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 230-239 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, /* 240-249 */
	0x00, 0x00, 0x00, 0x00, 0x00, 0x00}; /* 250-255 */

#define LF 0x0a
#define TAB 0x09
#define ZERO 0x30
#define ONE 0x31

int skipNtabs(unsigned char *x, int i, int n)
{
    int count = 0;
    while (count < n) {
	while (x[i] != TAB) i++;
	count++;
	i++;
    }
    return i;
}

int moveLeftNtabs(unsigned char *x, int i, int n)
{
    int count = 0;
    while (count < n) {
	while (x[i] != TAB) i--;
	count++;
	i--;
    }
    return ++i;
}

void getChromosome(unsigned char *x, int a, int b, char *chrom)
{
    int i, j;
    for (i = a, j = 0; i <= b; i++, j++) chrom[j] = x[i];
    chrom[j] = '\0';
}

/* from pegas */
static int raw2int(unsigned char *x, int a, int b)
{
	int i, k = 1, ans = 0;

	for (i = b; i >= a; i--, k *= 10)
		ans += ((int)x[i] - 48) * k;

	return ans;
}

SEXP vcf2DNAbinFromRaw(SEXP x, SEXP refgenome, SEXP sizeofheader, SEXP individual)
{
    unsigned char *xp, *seq, ALT, REF, base;
    int j, k, a, b, pos, Ninds = -8, eol;

    PROTECT(x = coerceVector(x, RAWSXP));
    PROTECT(refgenome = coerceVector(refgenome, VECSXP));
    PROTECT(sizeofheader = coerceVector(sizeofheader, INTSXP));
    PROTECT(individual = coerceVector(individual, INTSXP));

    SEXP res = duplicate(refgenome);
    //Rprintf("\nNAMED(refgenome) = %d\n", NAMED(refgenome));
    //Rprintf("NAMED(res) = %d\n", NAMED(res));

    int N = LENGTH(x);
    int NSEQ = LENGTH(refgenome); // nb of chromosomes
    SEXP seqnames = getAttrib(refgenome, R_NamesSymbol);
    int i = INTEGER(sizeofheader)[0];
    int ind = INTEGER(individual)[0];

    xp = RAW(x);
    char chrom[500];

    /* find number of individuals */
    for (k = i; xp[k] != LF; k++) if (xp[k] == TAB) Ninds++;
    if (Ninds < 1) error("apparently no genotypes in VCF file");
    if (Ninds < ind)
	error("argument 'individual' larger than the number of genotypes");

    for (;; i = eol + 1) {
	R_CheckUserInterrupt();
	if (i >= N) break;
	eol = i; // find end of line
	/* NOTE: the last locus line in the VCF file is not terminated
	   by a LF, so the checks below should be safe. */
	while (eol < N && xp[eol] != LF) eol++;
	if (eol > N) break;
	j = moveLeftNtabs(xp, eol, Ninds - ind + 1);
	/* check that the two alleles are not the REF allele */
	if (xp[j + 1] == ZERO && xp[j + 3] == ZERO) continue;
	/* check that the locus is an SNP */
	a = skipNtabs(xp, i, 3);
	if (xp[a + 1] != TAB) continue;
	if (xp[a + 3] != TAB) continue;
	REF = xp[a]; ALT = xp[a + 2];
	/* get the chromosome */
	b = i;
	while (xp[b] != TAB) b++;
	b--;
	getChromosome(xp, i, b, chrom);
	for (k = 0; k < NSEQ; k++)
	    if (strcmp(CHAR(STRING_ELT(seqnames, k)), chrom) == 0) break;
	if (k == NSEQ) error("mismatch between chromosome names");
	seq = RAW(VECTOR_ELT(res, k));
	/* get position of the locus */
	a = b += 2;
	while (xp[b] != TAB) b++;
	b--;
	pos = raw2int(xp, a, b);
	if (xp[j + 1] == ONE && xp[j + 3] == ONE) { // homozygote with ALT allele
	    seq[pos - 1] = tab_trans[ALT];
	} else { // heterozygote
	    base = tab_trans[REF] | tab_trans[ALT];
	    base &= 0xf7; /* set the 5th bit to 0 */
	    seq[pos - 1] = base;
	}
    }
    UNPROTECT(4);
    return res;
}

/* returns 8 if the base is known surely, 0 otherwise */
#define KnownBase(a) (a & 8)

/* returns 1 if the base is R, M, W, S, K, or Y; 0 otherwise */
#define isSNP(a) (a == 48 || a == 80 || a == 96 || a == 144 || a == 160 || a == 192)

SEXP seqBinning_DNAbin(SEXP x, SEXP binsize)
{
    int i, j, k, s, nseq, seql, rseql, next;
    unsigned char *xr, *rseq;
    SEXP obj, seq;
    PROTECT(x = coerceVector(x, VECSXP));
    nseq = LENGTH(x);
    PROTECT(binsize = coerceVector(binsize, INTSXP));
    s = INTEGER(binsize)[0];
    PROTECT(obj = allocVector(VECSXP, nseq));

    /* i: index of the sequences (in the lists)
       k: index of the base of the i-th sequence
       j: index of base of the output sequence */
    for (i = 0; i < nseq; i++) {
	xr = RAW(VECTOR_ELT(x, i));
	seql = LENGTH(VECTOR_ELT(x, i));
	rseql = seql/s;
	if (seql != rseql * s) rseql++;
	PROTECT(seq = allocVector(RAWSXP, rseql));
	rseq = RAW(seq);
	memset(rseq, 0x18, rseql * sizeof(unsigned char)); // initialize with T
	k = 0; next = s; j = 0;
	while (k < seql) {
	    if (! KnownBase(xr[k])) {
		if (isSNP(xr[k])) {
		    rseq[j] = 0x50; // K
		    k = next;
		    next += s;
		    j++;
		} else k++; // N, ...
	    } else k++; // A, C, G, or T
	    if (k == next) {
		next += s;
		j++;
	    }
	}
	SET_VECTOR_ELT(obj, i, seq);
	UNPROTECT(1);
    }

    UNPROTECT(3);
    return obj;
}
