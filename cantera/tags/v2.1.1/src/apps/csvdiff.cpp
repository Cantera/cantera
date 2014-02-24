/*
 *  csvdiff File1.csv File2.csv
 *
 *  Compares the variable values in two Excel formatted
 *  comma separated files.
 *  The comparison is done using a weighted norm basis.
 *
 *  The two files should be basically equal. However, File1.csv is
 *  taken as the reference file, that has precedence, when there is
 *  something to be decided upon.
 *
 *  Arguments:
 *   -h = prints this usage information
 *
 *  Shell Return Values
 *    0 = Comparison was successful
 *    1 = One or more nodal values failed the comparison
 *    2 = One or more of the header values failed the comparison
 *    3 = Apples to oranges, the files can not even be compared against
 *        one another.
 */

#include <math.h>
#include <float.h>
#include <algorithm>
#include <cstdlib>
#include "cantera/base/config.h"
#ifndef _MSC_VER
#include <unistd.h>
#else
#include <string>
#endif
using namespace std;

#if defined(__CYGWIN__)
#include <getopt.h>
#endif

#include "mdp_allo.h"
//#include "cantera/base/mdp_allo.h"
#include "tok_input_util.h"

int Debug_Flag = true;
double grtol = 1.0E-3;
double gatol = 1.0E-9;

#define RT_PASSED  0
#define RT_FAILED_COL 1
#define RT_FAILED_HDR 2
#define RT_FAILED_OTHER 3

/*
 * First iteration towards getting this variable
 */
int Max_Input_Str_Ln = MAX_INPUT_STR_LN;
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#ifdef _MSC_VER
/*
 * Windows doesn't have getopt(). This is an incomplete version that
 * does enough to handle required functionality.
 */
int optind = -1;
char* optarg = 0;

int getopt(int argc, char** argv, const char*)
{
    static int currArg = 1;
    static int currOptInd = 1;
    string tok;
    static int charPos = 0;
    int rc = -1;
    if (currArg >= argc) {
        optarg = 0;
        return -rc;
    }
    tok = string(argv[currArg]);
    currOptInd = currArg+1;
    if (currOptInd > argc - 1) {
        currOptInd = -1;
        optarg = 0;
    } else {
        optarg = argv[currArg+1];
    }
    size_t len = strlen(tok.c_str());
    if (charPos == 0) {
        bool found = false;
        do {
            tok = string(argv[currArg]);
            len = strlen(tok.c_str());
            if (len > 1 && tok[0] == '-') {
                found = true;
                charPos = 1;
                if (len > 2 && tok[1] == '-') {
                    charPos = 2;
                }
            } else {
                if (optind == -1) {
                    optind = currArg;
                }
            }
            if (!found) {
                if (currArg < (argc-1)) {
                    currArg++;
                } else {
                    optarg = 0;
                    return -1;
                }
            }
        } while (!found);
    }

    rc = tok[charPos];
    if (charPos < static_cast<int>(len - 1)) {
        charPos++;
    } else {
        charPos = 0;
    }
    return rc;
}

#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int diff_double(double d1, double d2, double rtol, double atol)

/*
 * Compares 2 doubles. If they are not within tolerance, then this
 * function returns true.
 */
{
    if (fabs(d1-d2) > (atol + rtol * 0.5 * (fabs(d1) + fabs(d2)))) {
        return 1;
    }
    return 0;
}

static int diff_double_slope(double d1, double d2, double rtol,
                             double atol, double xtol, double slope1, double slope2)

/*
 * Compares 2 doubles. If they are not within tolerance, then this
 * function returns true.
 */
{
    double atol2 = xtol*(fabs(slope1) + fabs(slope2));
    if (fabs(d1-d2) > (atol + atol2 + rtol * 0.5 * (fabs(d1) + fabs(d2)))) {
        return 1;
    }
    return 0;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static double calc_rdiff(double d1, double d2, double rtol, double atol)

/*
 * Calculates the relative difference using a fuzzy comparison
 */

{
    double rhs, lhs;
    rhs = fabs(d1-d2);
    lhs = atol + rtol * 0.5 * (fabs(d1) + fabs(d2));
    return rhs/lhs;
}

/*****************************************************************************/
/*
 * breakStrCommas():
 *    This routine will break a character string into stringlets according
 *    to the placement of commas. The commas are replaced by null
 *    characters.
 *
 *    Argument:
 *      str => original string. On exit, this string will have been altered.
 *      strlets -> Vector of pointers to char *.  The vector has a size
 *                 larger than or equal to maxPieces.
 *      maxPieces -> largest number of pieces to divide the string into.
 *
 *    Return:
 *      This returns the number of pieces that the string is actually
 *      broken up into.
 */

static int breakStrCommas(char* str, char** strlets, int maxPieces)
{
    int numbreaks = 0;
    if (strlets) {
        strlets[0] = str;
        if (str) {
            char* cptr = str;
            char* cetn = NULL;
            do {
                cetn = strchr(cptr, (int) ',');
                if (cetn) {
                    numbreaks++;
                    cptr = cetn + 1;
                    strlets[numbreaks] = cptr;
                    *cetn = '\0';
                }
            } while (cetn && (numbreaks < (maxPieces - 1)));
        }
    }
    return numbreaks + 1;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#define LT_NULLLINE  0
#define LT_TITLELINE 1
#define LT_COLTITLE  2
#define LT_DATALINE  3
/*
 * get_sizes()
 *
 * This routine obtains the sizes of the various elements of the file
 * by parsing the file.
 *   (HKM: Note, this file could use some work. However, it's always
 *         going to be heuristic)
 *
 *   Arguments:
 *
 *    fp           = File pointer
 *    nTitleLines = Number of title lines
 *    nColTitleLines = Number of column title lines
 *    nCol        = Number of columns -> basically equal to the
 *                  number of variables
 *    nDataRows   = Number of rows of data in the file
 *
 */

static void get_sizes(FILE* fp, int& nTitleLines, int& nColTitleLines,
                      int& nCol, int& nDataRows, int** ColIsFloat_ptr)
{
    int nScanLinesMAX = 100;
    int nScanLines = nScanLinesMAX;
    int retn, i, j;
    int maxCommas = 0;
    TOKEN fieldToken;
    char* scanLine = mdp_alloc_char_1(MAX_INPUT_STR_LN+1, '\0');
    int* numCommas = mdp_alloc_int_1(nScanLinesMAX, -1);
    int* ColIsFloat = *ColIsFloat_ptr;

    /*
     * Rewind the file
     */
    rewind(fp);
    /*
     * Read the scan lines
     */
    for (i = 0; i < nScanLinesMAX; i++) {
        retn = read_line(fp, scanLine, 0);
        if (retn == -1) {
            nScanLines = i;
            break;
        }
        /*
         * Strip a trailing comma from the scanline -
         *  -> These are not significant
         */
        int ccount = static_cast<int>(strlen(scanLine));
        if (ccount > 0) {
            if (scanLine[ccount-1] == ',') {
                scanLine[ccount-1] = '\0';
            }
        }
        /*
         * Count the number of commas in the line
         */
        char* cptr = scanLine;
        char* cetn = NULL;
        numCommas[i] = 0;
        do {
            cetn = strchr(cptr, (int) ',');
            if (cetn) {
                numCommas[i]++;
                cptr = cetn + 1;
            }
        } while (cetn);
        if (i > 1) {
            if (maxCommas < numCommas[i]) {
                maxCommas = numCommas[i];
            }
        }
    }
    /*
     * set a preliminary value of nCol
     */
    nCol = maxCommas + 1;
    if (nScanLines == 0) {
        nCol = 0;
    }
    char** strlets = (char**) mdp_alloc_ptr_1(maxCommas+1);

    /*
     * Figure out if each column is a text or float
     */
    rewind(fp);
    for (i = 0; i < nScanLines; i++) {
        retn = read_line(fp, scanLine, 0);
        int ncolsFound = breakStrCommas(scanLine, strlets, nCol);
        if (ncolsFound == (maxCommas + 1)) {
            for (j = 0; j < ncolsFound; j++) {
                char* fieldStr = strlets[j];
                fillTokStruct(&fieldToken, fieldStr);
                if (fieldToken.ntokes != 1) {
                    break;
                }
                bool rerr = false;
                (void) tok_to_double(&fieldToken, DBL_MAX,
                                     -DBL_MAX, 0.0, &rerr);
                if (!rerr) {
                    ColIsFloat[j] = true;
                }
            }

        }
    }



    int doingLineType = LT_TITLELINE;
    rewind(fp);
    for (i = 0; i < nScanLines; i++) {
        retn = read_line(fp, scanLine, 0);
        /*
         * Strip a trailing comma from the scanline -
         *  -> These are not significant
         */
        int ccount = static_cast<int>(strlen(scanLine));
        if (ccount > 0) {
            if (scanLine[ccount-1] == ',') {
                scanLine[ccount-1] = '\0';
            }
        }
        int ncolsFound = breakStrCommas(scanLine, strlets, nCol);

        if (doingLineType == LT_TITLELINE) {
            if (numCommas[i] == maxCommas) {
                doingLineType = LT_COLTITLE;
                nTitleLines = i;
            }
        }

        if (doingLineType == LT_COLTITLE) {
            bool goodDataLine = true;
            bool rerr = false;
            for (j = 0; j < ncolsFound; j++) {
                char* fieldStr = strlets[j];
                fillTokStruct(&fieldToken, fieldStr);
                if (fieldToken.ntokes != 1) {
                    goodDataLine = false;
                    break;
                }
                if ((ColIsFloat[j]) == 1) {
                    (void) tok_to_double(&fieldToken, DBL_MAX,
                                         -DBL_MAX, 0.0, &rerr);
                    if (rerr) {
                        goodDataLine = false;
                        break;
                    }
                }
            }
            if (goodDataLine) {
                doingLineType = LT_DATALINE;
            }
            nColTitleLines = i - nTitleLines;
        }
        if (doingLineType == LT_DATALINE) {
            break;
        }
    }


    /*
     * Count the total number of lines in the file
     */
    if (doingLineType == LT_DATALINE) {
        for (i = nColTitleLines + nTitleLines; ; i++) {
            retn = read_line(fp, scanLine, 0);
            if (retn == -1) {
                nDataRows = i - nColTitleLines - nTitleLines + 1;
                break;
            }
            /*
             * Strip a trailing comma from the scanline -
             *  -> These are not significant
             */
            int ccount = static_cast<int>(strlen(scanLine));
            if (ccount > 0) {
                if (scanLine[ccount-1] == ',') {
                    scanLine[ccount-1] = '\0';
                }
            }
            int ncolsFound = breakStrCommas(scanLine, strlets, nCol);
            bool goodDataLine = true;
            bool rerr = false;
            for (j = 0; j < ncolsFound; j++) {
                char* fieldStr = strlets[j];
                fillTokStruct(&fieldToken, fieldStr);
                if (fieldToken.ntokes != 1) {
                    goodDataLine = false;
                    break;
                }
                if (ColIsFloat[j] == 1) {
                    (void) tok_to_double(&fieldToken, DBL_MAX,
                                         -DBL_MAX, 0.0, &rerr);
                    if (rerr) {
                        goodDataLine = false;
                        break;
                    }
                }
            }
            if (! goodDataLine) {
                doingLineType = LT_NULLLINE;
                nDataRows = i - nColTitleLines - nTitleLines + 1;
                break;
            }
        }
    }
    mdp_safe_free((void**) &strlets);
    mdp_safe_free((void**) &scanLine);
    mdp_safe_free((void**) &numCommas);
    return;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void
read_title(FILE* fp, char** *title, int nTitleLines)
{
    int retn;
    *title = (char**) mdp_alloc_ptr_1(nTitleLines);
    char* scanLine = mdp_alloc_char_1(Max_Input_Str_Ln + 1, '\0');
    for (int i = 0; i < nTitleLines ; i++) {
        retn = read_line(fp, scanLine, 0);
        if (retn >= 0) {
            /*
             * Strip a trailing comma from the scanline -
             *  -> These are not significant
             */
            int ccount = static_cast<int>(strlen(scanLine));
            if (ccount > 0) {
                if (scanLine[ccount-1] == ',') {
                    scanLine[ccount-1] = '\0';
                }
            }
            (*title)[i] = mdp_copy_string(scanLine);
        }
    }
    mdp_safe_free((void**) &scanLine);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void
read_colTitle(FILE* fp,  char**** ColMLNames_ptr, int  nColTitleLines, int nCol)
{
    int retn, j;
    *ColMLNames_ptr = (char***) mdp_alloc_ptr_1(nCol);
    char** *ColMLNames = *ColMLNames_ptr;
    char* scanLine = mdp_alloc_char_1(Max_Input_Str_Ln + 1, '\0');
    char** strlets = (char**) mdp_alloc_ptr_1(nCol+1);
    if (nColTitleLines > 0) {
        for (int i = 0; i < nColTitleLines ; i++) {
            retn = read_line(fp, scanLine, 0);
            if (retn >= 0) {
                /*
                 * Strip a trailing comma from the scanline -
                 *  -> These are not significant
                 */
                int ccount = static_cast<int>(strlen(scanLine));
                if (ccount > 0) {
                    if (scanLine[ccount-1] == ',') {
                        scanLine[ccount-1] = '\0';
                    }
                }
                int ncolsFound = breakStrCommas(scanLine, strlets, nCol);
                ColMLNames[i] = mdp_alloc_VecFixedStrings(nCol, MAX_TOKEN_STR_LN+1);
                for (j = 0; j < ncolsFound; j++) {
                    strip(strlets[j]);
                    strcpy(ColMLNames[i][j], strlets[j]);
                }
            }
        }
    } else {
        ColMLNames[0] = mdp_alloc_VecFixedStrings(nCol, MAX_TOKEN_STR_LN+1);
        for (j = 0; j < nCol; j++) {
            char cbuff[256];
            sprintf(cbuff, "Col_%d", j+1);
            strcpy(ColMLNames[0][j], cbuff);
        }
    }
    mdp_safe_free((void**) &scanLine);
    mdp_safe_free((void**) &strlets);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static double get_atol(const double* values, const int nvals,
                       const double atol)
{
    int i;
    double sum = 0.0, retn;
    if (nvals <= 0) {
        return gatol;
    }
    for (i = 0; i < nvals; i++) {
        retn = values[i];
        sum += retn * retn;
    }
    sum /= nvals;
    retn = sqrt(sum);
    return (retn + 1.0) * atol;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void
read_values(FILE* fp, double** NVValues, char** *NSValues, int nCol, int nDataRows,
            int* ColIsFloat)
{
    char** strlets = (char**) mdp_alloc_ptr_1(nCol+1);
    char* scanLine = mdp_alloc_char_1(Max_Input_Str_Ln + 1, '\0');
    TOKEN fieldToken;
    double value;
    int retn, j;
    for (int i = 0; i < nDataRows; i++) {
        retn = read_line(fp, scanLine, 0);
        if (retn == -1) {
            break;
        }
        /*
         * Strip a trailing comma from the scanline -
         *  -> These are not significant
         */
        int ccount = static_cast<int>(strlen(scanLine));
        if (ccount > 0) {
            if (scanLine[ccount-1] == ',') {
                scanLine[ccount-1] = '\0';
            }
        }
        int ncolsFound = breakStrCommas(scanLine, strlets, nCol);
        bool goodDataLine = true;
        bool rerr = false;
        for (j = 0; j < ncolsFound; j++) {
            char* fieldStr = strlets[j];
            NSValues[j][i] = mdp_copy_string(strlets[j]);
            fillTokStruct(&fieldToken, fieldStr);
            if (fieldToken.ntokes != 1) {
                goodDataLine = false;
                break;
            }
            if (ColIsFloat[j]) {
                value = tok_to_double(&fieldToken, DBL_MAX,
                                      -DBL_MAX, 0.0, &rerr);
                if (rerr) {
                    goodDataLine = false;
                    break;
                }
                NVValues[j][i] = value;
            }
        }
        if (! goodDataLine) {
            break;
        }
    }
    mdp_safe_free((void**) &strlets);
    mdp_safe_free((void**) &scanLine);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void print_usage()
{
    printf("\t\n");
    printf("  csvdiff [-h] [-a atol] [-r rtol] File1.csv File2.csv\n");
    printf("\t\n");
    printf("\tCompares the variable values in two Excel formatted "
           "comma separated files.\n");
    printf("\tThe comparison is done using a weighted norm basis.\n");
    printf("\t\n");
    printf("\tThe two files should be basically equal. However, File1.csv is\n");
    printf("\ttaken as the reference file that has precedence, when there is\n");
    printf("\tsomething to be decided upon.\n");
    printf("\t\n");
    printf("\t Arguments:\n");
    printf("\t  -h      = Usage info\n");
    printf("\t  -a atol = Set absolute tolerance parameter - default = 1.0E-9\n");
    printf("\t  -r rtol = Set relative tolerance parameter - default = 1.0E-3\n");
    printf("\t\n");
    printf("\t Shell Return Values:\n");
    printf("\t   0 = Comparison was successful\n");
    printf("\t   1 = One or more nodal values failed the comparison\n");
    printf("\t   2 = One or more header values failed the comparison\n");
    printf("\t   3 = Apples to oranges, the files can not even be compared against\n");
    printf("\t       one another.\n");
    printf("\t\n");
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int main(int argc, char* argv[])

/*
 *                        main driver for csvdiff.
 */
{
    int opt_let;
    char*  fileName1=NULL, *fileName2=NULL; /* Names of the csv files */
    FILE*  fp1=NULL, *fp2=NULL;
    int    nTitleLines1 = 0, nTitleLines2 = 0;
    int    nColTitleLines1 = 0, nColTitleLines2 = 0;
    int    nCol1 = 0,  nCol2 = 0, nColMAX = 0, nColcomparisons = 0;
    int    nDataRows1 = 0, nDataRows2 = 0;
    char**  title1 = 0,  **title2 = 0;
    int**    compColList = NULL;
    char**   *ColMLNames1 = NULL, *** ColMLNames2 = NULL;
    char**   ColNames1 = NULL, **ColNames2 = NULL;
    double** NVValues1 = NULL, **NVValues2 = NULL;
    char**   *NSValues1 = NULL, *** NSValues2 = NULL;
    int*    ColIsFloat1 = NULL, *ColIsFloat2 = NULL;
    double* curVarValues1 = NULL, *curVarValues2 = NULL;
    char** curStringValues1 = NULL, **curStringValues2 = NULL;
    int    i, j, ndiff, jmax=0, i1, i2, k;
    bool found;
    double max_diff, rel_diff;
    int    testPassed = RT_PASSED;
    double atol_j, atol_arg = 0.0, rtol_arg = 0.0;

    /********************** BEGIN EXECUTION ************************************/
    int id = 0;
    int id2 = 0;
    char* ggg = 0;
    char* rrr = 0;
    /*
    * Interpret command line arguments
    */
    /* Loop over each command line option */
    while ((opt_let = getopt(argc, argv, "ha:r:")) != EOF) {

        /* case over the option letter */
        switch (opt_let) {

        case 'h':
            /* Usage info was requested */
            print_usage();
            exit(0);

        case 'a':
            /* atol parameter */

            ggg = optarg;
            //printf("a = %s\n", ggg);
            id = sscanf(ggg,"%lg", &atol_arg);
            if (id != 1) {
                printf(" atol param bad: %s\n", ggg);
                exit(-1);
            }
            gatol = atol_arg;
            break;

        case 'r':
            /* rtol parameter */

            rrr = optarg;
            //printf("r = %s\n", ggg);
            id2 = sscanf(rrr,"%lg", &rtol_arg);
            if (id2 != 1) {
                printf(" rtol param bad: %s\n", rrr);
                exit(-1);
            }
            grtol = rtol_arg;
            break;


        default:
            /* Default case. Error on unknown argument. */
            printf("default called opt_let = %c\n", opt_let);
            fprintf(stderr, "ERROR in command line usuage:\n");
            print_usage();
            return 0;
        } /* End "switch(opt_let)" */

    } /* End "while((opt_let=getopt(argc, argv, "i")) != EOF)" */

    if (optind !=  argc-2) {
        print_usage();
        exit(-1);
    } else {
        fileName1 = argv[argc-2];
        fileName2 = argv[argc-1];
    }

    /*
     *      Print Out Header
     */
    printf("\n");
    printf("----------------------------------------------------------\n");
    printf("csvdiff: CSVFile comparison utility program\n");
    printf("         Harry K. Moffat Div. 9114 Sandia National Labs\n");
    printf("         \n");
    printf("         First  CSV File = %s\n", fileName1);
    printf("         Second CSV file = %s\n", fileName2);
    printf("\n");
    printf("         Absolute tol = %g\n", gatol);
    printf("         Relative tol = %g\n", grtol);
    printf("----------------------------------------------------------\n");
    printf("\n");

    /*
    *  Open up the two ascii Files #1 and #2
    */
    if (!(fp1 = fopen(fileName1, "r"))) {
        fprintf(stderr,"Error opening up file1, %s\n", fileName1);
        exit(-1);
    }
    if (!(fp2 = fopen(fileName2, "r"))) {
        fprintf(stderr, "Error opening up file2, %s\n", fileName2);
        exit(-1);
    }

    ColIsFloat1 = mdp_alloc_int_1(200, 0);
    ColIsFloat2 = mdp_alloc_int_1(200, 0);
    /*
     *   Obtain the size of the problem information: Compare between files.
     */

    get_sizes(fp1, nTitleLines1, nColTitleLines1, nCol1, nDataRows1, &ColIsFloat1);
    if (nCol1 == 0) {
        printf("Number of columns in file %s is zero\n", fileName1);
        testPassed = RT_FAILED_OTHER;
        exit(RT_FAILED_OTHER);
    }
    if (nDataRows1 == 0) {
        printf("Number of data rows in file %s is zero\n", fileName1);
        testPassed = RT_FAILED_OTHER;
        exit(RT_FAILED_OTHER);
    }

    get_sizes(fp2, nTitleLines2, nColTitleLines2, nCol2, nDataRows2, &ColIsFloat2);
    if (nCol2 == 0) {
        printf("Number of columns in file %s is zero\n", fileName2);
        testPassed = RT_FAILED_OTHER;
        exit(RT_FAILED_OTHER);
    }
    if (nDataRows2 == 0) {
        printf("Number of data rows in file %s is zero\n", fileName2);
        testPassed = RT_FAILED_OTHER;
        exit(RT_FAILED_OTHER);
    }

    if (nTitleLines1 != nTitleLines2) {
        printf("Number of Title Lines differ:, %d %d\n",nTitleLines1, nTitleLines2);
        testPassed = RT_FAILED_OTHER;
    } else if (Debug_Flag) {
        printf("Number of Title Lines in each file = %d\n", nTitleLines1);
    }
    if (nColTitleLines1 != nColTitleLines2) {
        printf("Number of Column title lines differ:, %d %d\n", nColTitleLines1,
               nColTitleLines2);
        testPassed = RT_FAILED_OTHER;
    } else if (Debug_Flag) {
        printf("Number of column title lines in each file = %d\n", nColTitleLines1);
    }

    /*
     * Right now, if the number of data rows differ, we will punt.
     * Maybe later we can do something more significant
     */
    int nDataRowsMIN = min(nDataRows1, nDataRows2);
    int nDataRowsMAX = max(nDataRows1, nDataRows2);
    if (nDataRows1 != nDataRows2) {
        printf("Number of Data rows in file1, %d, is different than file2, %d\n",
               nDataRows1, nDataRows2);
    } else {
        printf("Number of Data rows in both files = %d\n", nDataRowsMIN);
    }

    rewind(fp1);
    rewind(fp2);
    read_title(fp1, &title1, nTitleLines1);
    read_title(fp2, &title2, nTitleLines2);

    if (nTitleLines1 > 0 || nTitleLines2 > 0) {
        int n = min(nTitleLines1, nTitleLines2);
        for (i = 0; i < n; i++) {
            if (strcmp(title1[i], title2[i]) != 0) {
                printf("Title Line %d differ:\n\t\"%s\"\n\t\"%s\"\n", i, title1[i], title2[i]);
                testPassed = RT_FAILED_HDR;
            } else if (Debug_Flag) {
                printf("Title Line %d for each file: \"%s\"\n", i, title1[i]);
            }
        }
        if (nTitleLines1 != nTitleLines2) {
            printf("Number of Title Lines differ: %d %d\n", nTitleLines1, nTitleLines2);
            testPassed = RT_FAILED_HDR;
        }
    } else {
        if (nTitleLines1 != nTitleLines2) {
            if (nTitleLines1) {
                printf("Titles differ: title for first file: \"%s\"\n",
                       title1[0]);
                testPassed = RT_FAILED_HDR;
            }
            if (nTitleLines2) {
                printf("Titles differ: title for second file: \"%s\"\n",
                       title2[0]);
            }
            testPassed = RT_FAILED_HDR;
        }
    }

    /*
     *  Get the number of column variables in each file
     */

    if (nCol1 != nCol2) {
        printf("Number of column variables differ:, %d %d\n",
               nCol1, nCol2);
        testPassed = RT_FAILED_OTHER;
    } else if (Debug_Flag) {
        printf("Number of column variables in both files = %d\n",
               nCol1);
    }

    /*
     * Read the names of the column variables
     */
    read_colTitle(fp1, &ColMLNames1,  nColTitleLines1, nCol1);
    read_colTitle(fp2, &ColMLNames2,  nColTitleLines2, nCol2);
    ColNames1 = ColMLNames1[0];
    ColNames2 = ColMLNames2[0];

    /*
     * Do a Comparison of the names to find the maximum number
     * of matches.
     */
    nColMAX = max(nCol1, nCol2);

    compColList = mdp_alloc_int_2(nColMAX, 2, -1);
    nColcomparisons = 0;
    for (i = 0; i < nCol1; i++) {
        found = false;
        for (j = 0; j < nCol2; j++) {
            if (!strcmp(ColNames1[i], ColNames2[j])) {
                compColList[nColcomparisons][0] = i;
                compColList[nColcomparisons][1] = j;
                nColcomparisons++;
                found = true;
                break;
            }
        }
        if (!found) {
            printf("csvdiff WARNING Variable %s (%d) in first file not found"
                   " in second file\n", ColNames1[i], i);
            testPassed = RT_FAILED_OTHER;
        }
    }
    for (j = 0; j < nCol2; j++) {
        found = false;
        for (i = 0; i <  nColcomparisons; i++) {
            if (compColList[i][1] == j) {
                found = true;
            }
        }
        if (! found) {
            printf("csvdiff WARNING Variable %s (%d) in second file "
                   "not found in first file\n",
                   ColNames2[j], j);
            testPassed = RT_FAILED_OTHER;
        }
    }

    /*
     *  Allocate storage for the column variables
     */
    NVValues1 = mdp_alloc_dbl_2(nCol1, nDataRowsMAX, 0.0);
    NVValues2 = mdp_alloc_dbl_2(nCol2, nDataRowsMAX, 0.0);

    /*
     *  Allocate storage for the column variables
     */
    NSValues1 = (char***) mdp_alloc_ptr_2(nCol1, nDataRowsMAX);
    NSValues2 = (char***) mdp_alloc_ptr_2(nCol2, nDataRowsMAX);

    /*
     *  Read in the values to the arrays
     */
    read_values(fp1, NVValues1, NSValues1, nCol1, nDataRows1, ColIsFloat1);
    read_values(fp2, NVValues2, NSValues2, nCol2, nDataRows2, ColIsFloat2);

    /*
     * Compare the solutions in each file
     */
    int method = 1;
    double slope1, slope2, xatol;
    int notOK;
    for (k = 0; k < nColcomparisons; k++) {

        i1 =  compColList[k][0];
        i2 =  compColList[k][1];
        bool doFltComparison = true;
        if (!ColIsFloat1[i1]) {
            doFltComparison = false;
            jmax = -1;
        }
        if (!ColIsFloat2[i2]) {
            doFltComparison = false;
            jmax = -1;
        }
        curStringValues1 = NSValues1[i1];
        curStringValues2 = NSValues2[i2];
        max_diff = 0.0;
        ndiff = 0;
        if (doFltComparison) {
            curVarValues1 = NVValues1[i1];
            curVarValues2 = NVValues2[i2];
            atol_j =             get_atol(curVarValues1, nDataRows1, gatol);
            atol_j = min(atol_j, get_atol(curVarValues2, nDataRows2, gatol));
            for (j = 0; j < nDataRowsMIN; j++) {

                slope1 = 0.0;
                slope2 = 0.0;
                xatol = fabs(grtol * (NVValues1[0][j] - NVValues1[0][j-1]));
                if (j > 0 && k > 0) {
                    slope1 = (curVarValues1[j] - curVarValues1[j-1])/
                             (NVValues1[0][j] - NVValues1[0][j-1]);
                    slope2 = (curVarValues2[j] - curVarValues2[j-1])/
                             (NVValues2[0][j] - NVValues2[0][j-1]);
                }
                if (method) {
                    notOK = diff_double_slope(curVarValues1[j], curVarValues2[j],
                                              grtol, atol_j, xatol, slope1, slope2);
                } else {
                    notOK = diff_double(curVarValues1[j], curVarValues2[j],
                                        grtol, atol_j);
                }
                if (notOK) {
                    ndiff++;
                    rel_diff = calc_rdiff((double) curVarValues1[j],
                                          (double) curVarValues2[j], grtol, atol_j);
                    if (rel_diff > max_diff) {
                        jmax = j;
                        max_diff = rel_diff;
                    }
                    if (ndiff < 10) {
                        printf("\tColumn variable %s at data row %d ", ColNames1[i1], j + 1);
                        printf(" differ: %g %g\n", curVarValues1[j],
                               curVarValues2[j]);
                    }
                }
            }
        } else {
            for (j = 0; j < nDataRowsMIN; j++) {
                strip(curStringValues1[j]);
                strip(curStringValues2[j]);
                notOK = false;
                if (strcmp(curStringValues1[j], curStringValues2[j])) {
                    notOK = true;
                    ndiff++;
                    if (ndiff < 10) {
                        printf("\tColumn String variable %s at data row %d ", ColNames1[i1], j + 1);
                        printf(" differ: %s %s\n", curStringValues1[j],
                               curStringValues2[j]);
                    }
                }
            }
        }

        if (nDataRowsMIN != nDataRowsMAX) {
            ndiff +=  nDataRowsMAX - nDataRowsMIN;
            if (ndiff < 10) {
                if (nDataRows1 > nDataRows2) {
                    for (j = nDataRowsMIN; j < nDataRowsMAX; j++) {
                        printf("\tColumn variable %s at data row %d ", ColNames1[i1], j + 1);
                        printf(" differ: %g      NA\n", curVarValues1[j]);
                    }
                } else {
                    for (j = nDataRowsMIN; j < nDataRowsMAX; j++) {
                        printf("\tColumn variable %s at data row %d ", ColNames1[i1], j + 1);
                        printf(" differ: NA     %g \n", curVarValues2[j]);
                    }
                }
            }
        }

        /*
         *  Print out final results of nodal variable test
         */

        if (ndiff > 0) {
            printf(
                "Column variable %s failed comparison test for %d occurrences\n",
                ColNames1[i1], ndiff);
            if (jmax >= 0) {
                printf("  Largest difference was at data row %d ", jmax + 1);
                printf(": %g %g\n", curVarValues1[jmax],  curVarValues2[jmax]);
            }
            testPassed = RT_FAILED_COL;
        } else if (Debug_Flag) {
            printf("Column variable %s passed\n",  ColNames1[i1]);
        }

    }

    return testPassed;

} /************END of main() *************************************************/
/*****************************************************************************/
