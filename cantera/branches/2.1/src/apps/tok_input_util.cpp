#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <ctype.h>

#include "tok_input_util.h"

static char DEFAULT_STR[8] = "default";
static char DELIMITERS[8] = " \t\n\f\r\v"; /* Defn of white space in isspace()
                                              - Used in tokenizing lines */
static char COM_CHAR = '!';   /* This is used as a comment character */
static char COM_CHAR2 = '#';  /* This is used as a 2nd comment character */
static char KEY_CHAR = '=';   /* This is used to separate the key_string
                                 from the rest of the line */
static int PrintInputFile = true;     /* Used to turn on and off the
                                      printing of the input file */

/*************** R O U T I N E S   I N   T H E   F I L E *******************
*
*    NAME                       TYPE            CALLED_BY
*--------------------------------------------------------------------
*    get_next_keyLine           bool            extern
*    tok_to_int                 int             extern
*    str_to_int                 int             extern, tok_to_int
*    tok_to_double              double          extern
*    str_to_double              double          extern,tok_to_double
*    tok_to_boolean             bool            extern
*    str_to_boolean             bool            extern,tok_to_boolean
*    tok_to_string              char *          extern
*
*    scan_for_int               int             extern
*    scan_for_double            double          extern
*    scan_for_string            char *          extern
*    scan_for_boolean           bool            extern
*    scan_for_line              int             extern
*    read_line                  int             scan_for_line,
*                                               get_next_keyLine
*    interpret_int              static bool     str_to_int + others
*    interpret_boolean          static bool     str_to_boolean
*    interpret_double           static bool     str_to_double
*    strip                      int             read_input_file,
*                                               look_for,
*                                               get_next_keyLine
*    read_string                static void     scan_for_line
*    stokenize                  int             fillTokStruct
*    outofbnds                  static bool     all
*    strmatch                   bool            extern, toktokmatch
*    strstrmatch                bool            extern
*    strtokmatch                bool            extern
*    toktokmatch                bool            extern, strtokmatch
*                                               strstrmatch
*    fillTokStruct              void            extern, strtokmatch
*                                               strstrmatch,
*                                               get_next_keyLine
*    copyTokStruct              void            extern
*
******************************************************************************/
/*
* Definitions of static functions:
*/

static bool outofbnds(const double, const double, const double);
static bool interpret_boolean(const char*, int*, const int);
static bool interpret_int(const char*, int*, const int, const int,
                          const int);
static bool interpret_double(const char*, double*, const double,
                             const double, const double);


/************ Member Functions for the TOKEN Structure ********************/

TOKEN::TOKEN(void) :
    orig_str(0),
    tok_str(0),
    ntokes(0)
{
    orig_str =  copy_string("");
    tok_str  =  copy_string("");
    tok_ptr[0] = orig_str;
}


TOKEN::TOKEN(const char* str) :
    orig_str(0),
    tok_str(0),
    ntokes(0)
{
    if (str == NULL) {
        orig_str = copy_string("");
        tok_str  = copy_string("");
        tok_ptr[0] = orig_str;
        ntokes = 0;
    } else {
        orig_str = copy_string(str);
        tok_str  = copy_string(str);
        ntokes = stokenize(tok_str, DELIMITERS, tok_ptr, MAXTOKENS);
    }
}

TOKEN::~TOKEN()
{
    if (orig_str) {
        free(orig_str);
    }
    orig_str = NULL;
    if (tok_str) {
        free(tok_str);
    }
    tok_str = NULL;
}

/**************************************************************************/

bool get_next_keyLine(FILE* ifp, TOKEN* keyLineTok, TOKEN* keyArgTok)
/*
 *   This routine reads the input file to obtain the next line of
 *   uncommented
 *   data. The results are returned in two TOKEN structures. keyLineTok
 *   contains the key Line (everything before the first equals sign).
 *   keyArgTok contains everything after the equals sign.
 *   Note - Either keyLineTok or keyArgTok may be the null token
 *                 (but not both)
 *
 *   The definition of a token structure, given in .h file,
 *   is as follows:
 *
 *      struct TOKEN {
 *        char  orig_str[MAX_INPUT_STR_LN + 1];
 *        char  tok_str[MAX_INPUT_STR_LN + 1];
 *        char *tok_ptr[MAXTOKENS];
 *        int   ntokes;
 *      };
 *   mdp_allo.h
 *     orig_str Contains the original string, unmodified.
 *     tok_str  Contains a modified version of the string,
 *              whose positions
 *              are pointed to by tok_ptr[i] values. It is usually not
 *              referenced directly.
 *     tok_ptr[i] Contains the i_th token of the original string. This
 *              is a stripped character string. 0 <=i <= i-1
 *     ntokes   Number of tokens in the string.
 *
 *
 *   Comments are denoted by either '!' or '#'.
 *   Everything after the comment character on a
 *   line is stripped first. The comment character can occur
 *   anywhere in the line.
 *
 *   Arguments to the keyLine are denoted by everything after a
 *   '=' character on the line.
 *
 *   Example:
 *   ---------------------------------------------------------------
 *   ! Jack and Jill went up the hill to fetch a pale of water
 *   ! Jack by nimble, Jack be swift; Jack jump over the candle stick
 *
 *   The meaning of life is =    36.243     24  136 Not even ! close
 *   -----------------------------------------------------------------
 *
 *   Then, the routine would return (amongst other things):
 *     keyLineTok->orig_str = "The meaning of life is"
 *     keyArgTok->orig_str = "36.243     24   36 Not even"
 *     keyArgTok->ntokes = 5
 *     keyArgTok->tok_ptr[0] = "36.243"
 *     keyArgTok->tok_ptr[1] = "24"
 *     keyArgTok->tok_ptr[2] = "136"
 *     keyArgTok->tok_ptr[3] = "Not"
 *     keyArgTok->tok_ptr[4] = "Even"
 *
 *   The function returns true if there is a next line to process.
 *   It returns false if an EOF is encountered.
 */
{
    int       retn_value, i;
    char      save_input[MAX_INPUT_STR_LN + 1];
    char*     token_start = NULL;

    /*
     *  Check the arguments to the routine. This routine needs to be
     *  supplied with valid pointers to files and spaces for storage
     *  of its output in the return tokens.
     */
    if (ifp == NULL || keyLineTok == NULL || keyArgTok == NULL) {
        fprintf(stderr, "get_next_keyLine ERROR, arguments are bad\n");
        return false;
    }

    /*
     * Read a chunk of text, either up to a newline from the file pointer,
     * ifp. If an EOF occurs, return without changing the input structure
     */
do_it_again:
    do {
        if ((retn_value = read_string(ifp, save_input, '\n')) < 0) {
            return false;
        }
        if (PrintInputFile)  {
            if (retn_value <=0) {
                printf("%s\n", save_input);
            } else {
                printf("%s",   save_input);
            }
        }
        for (i = 0; i < (int) strlen(save_input); i++) {
            if (save_input[i] == COM_CHAR || save_input[i] == COM_CHAR2) {
                save_input[i] = '\0';
                break;
            }
        }
    } while (strip(save_input) == 0);

    /*
    *  Discover whether there are arguments in the line
    *  and then separate the line into two
    */

    for (i = 0; i < (int) strlen(save_input); i++) {
        if (save_input[i] == KEY_CHAR) {
            save_input[i] = '\0';
            token_start = save_input + i + 1;
            break;
        }
    }

    /*
     *  Strip the two strings of leading and trailing white space.
     *  If both strings are now the null string (because the line
     *  consisted of a single '=' character for example),
     *  go back and get a new line.
     */

    i = strip(token_start);
    if (!strip(save_input)) if (!i) {
            goto do_it_again;
        }

    /*
     *  Now that we have two strings representing the Key String and
     *  associated arguments, process them into TOKEN structures.
     *  Note - if token_start still points to NULL, then fillTokStruct
     *         will fill keyArgTok with a "null token".
     */

    fillTokStruct(keyLineTok, save_input);
    fillTokStruct(keyArgTok,  token_start);
    return (true);
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

int tok_to_int(const TOKEN* tokPtr, const int maxVal, const int minVal,
               const int defaultVal, bool* error)
/*
 *      Interprets the first string of a TOKEN structure as an int.
 *      Returns the interpreted value as the return value.
 *      AnErrors condition is created if more than one token is found
 *      in the struct TOKEN.
 *      Bounds checking is done on the value before returning.  Value
 *      must be between the maxVal and minVal; it can equal the max or min
 *      value.
 *
 *      Certain ascii strings are checked for first (case is insignificant):
 *
 *                 String              Retn_Value (defined in <limits.h>
 *                 ---------        --------------
 *                  INT_MAX, max, all     INT_MAX
 *                  INT_MIN               INT_MIN
 *                  N/A, Not Available    INT_MIN
 *                  default or ""         defaultVal
 *
 *      A default may be specified on the command line. The absence of a
 *      default may also be specified by setting default_value to
 *      NO_DEFAULT_INT.
 *
 *      If there is an error, *error is set to true. *error isn't touched
 *      if there isn't an error.
 */
{
    if (tokPtr->ntokes == 0) {
        return str_to_int(DEFAULT_STR, maxVal, minVal, defaultVal, error);
    } else if (tokPtr->ntokes > 1) {
        (void) fprintf(stderr, "ERROR: tok_to_int, ntokes > 1: %s\n",
                       tokPtr->orig_str);
        *error = true;
    }
    return str_to_int(tokPtr->tok_ptr[0], maxVal, minVal, defaultVal, error);
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

int str_to_int(const char* int_string, const int maxVal, const int minVal,
               const int defaultVal, bool* error)
/*
 *      Interprets a stripped character string as an integer.
 *      Bounds checking is done on the value before returning.  Value
 *      must be between the max and min; it can equal the max or min value.
 *
 *      Certain ascii strings are checked for first (case is insignificant):
 *
 *                 String              Retn_Value (defined in <limits.h>
 *                 ---------        --------------
 *                  INT_MAX, max, all     INT_MAX
 *                  INT_MIN               INT_MIN
 *                  N/A, Not Available    INT_MIN
 *                  default               default_value
 *
 *      A default may be specified on the command line. The absence of a
 *      default may also be specified by setting default_value to
 *      NO_DEFAULT_INT.
 *
 *      If there is an error, *error is set to true. *error isn't touched
 *      if there isn't an error.
 */
{
    int retn_value, check = false;
    double LmaxVal, LminVal;
    if (defaultVal == NO_DEFAULT_INT) {
        check = true;
    }
    if (interpret_int(int_string, &retn_value, maxVal, minVal, defaultVal)) {
        if (check) {
            if (retn_value == NO_DEFAULT_INT) {
                (void) fprintf(stderr,
                               "ERROR: str_to_int: Default not allowed\n");
                *error = true;
            }
        }
        if (maxVal < INT_MAX) {
            LmaxVal = (double) maxVal + 0.01;
        } else {
            LmaxVal = (double) maxVal;
        }
        if (minVal > INT_MIN) {
            LminVal = (double) minVal - 0.01;
        } else {
            LminVal = (double) minVal;
        }
        if (outofbnds((double) retn_value, LmaxVal, LminVal)) {
            fprintf(stderr,
                    "ERROR: str_to_int outofbnds:\n\t\"%s\"\n",
                    int_string);
            fprintf(stderr,"\tmax = %d, min = %d\n", maxVal, minVal);
            *error = true;
        }
    } else {
        *error = true;
    }
    return retn_value;
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

double tok_to_double(const TOKEN* tokPtr, const double maxVal,
                     const double minVal, const double defaultVal,
                     bool* error)
/*
 *      Interprets the first string of a TOKEN structure as a double.
 *      Returns the interpreted value as the return value.
 *      Errors conditions are created if more than one token is found.
 *
 *      Bounds checking is then done on the value before returning.  Value
 *      must be between the max and min; it can equal the max or min.
 *
 *      Useful values for bounds (specified in <limits.h>:
 *          DBL_MAX = largest legitimate value of a double ~ 2.0E-308
 *         -DBL_MAX = smallest legitimate value of a double ~ 2.0E-308
 *          DBL_EPSILON = smallest value of a double that can be added to one
 *                        and produce a different number. ~ 2.0E-16
 *          DBL_MIN = smallest value of a double ~ 2.0E-308
 *      For example:
 *        If 0.0 is not a legitimate number for value, set min = DBL_MIN
 *        If value>=0.0 is legitimate, set min = 0.0
 *        if value<=100., set max = 100.
 *        If no range checking is required, set max = DBL_MAX, min = -DBL_MAX
 *
 *      Certain ascii strings are checked for first (case is insignificant):
 *
 *                 Token_String      Retn_Value ( specified in <limits.h>
 *                 ---------        ---------------------------------------
 *                  FLT_MAX, all          FLT_MAX
 *                  DBL_MAX, max          DBL_MAX
 *                  N/A                  -DBL_MAX
 *                  default or ""         default_value
 *                  small, DBL_MIN        DBL_MIN
 *                  DBL_EPSILON           DBL_EPSILON
 *
 *      A default may be specified by either specifying the value as
 *      "default" or by the absense of any tokens in the TOKEN struct.
 *      The absence of a default may also be specified by setting the
 *      value of default_value to NO_DEFAULT_DOUBLE.
 *
 *      If there is an error, *error is set to true. *error isn't touched
 *      if there isn't an error.
 */
{
    if (tokPtr->ntokes == 0) {
        return str_to_double(DEFAULT_STR, maxVal, minVal, defaultVal, error);
    } else if (tokPtr->ntokes > 1) {
        (void) fprintf(stderr, "ERROR: tok_to_double, ntokes > 1: %s\n",
                       tokPtr->orig_str);
        *error = true;
    }
    return str_to_double(tokPtr->tok_ptr[0], maxVal, minVal, defaultVal, error);
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

double str_to_double(const char* dbl_string, const double maxVal,
                     const double minVal, const double defaultVal,
                     bool* error)
/*
 *      Interprets a stripped character string as a double. Returns the
 *      interpreted value as the return value.
 *
 *      Bounds checking is then done on the value before returning.  Value
 *      must be between the max and min; it can equal the max or min.
 *
 *      Useful values for bounds (specified in <limits.h>:
 *          DBL_MAX = largest legitimate value of a double ~ 2.0E-308
 *         -DBL_MAX = smallest legitimate value of a double ~ 2.0E-308
 *          DBL_EPSILON = smallest value of a double that can be added to one
 *                        and produce a different number. ~ 2.0E-16
 *          DBL_MIN = smallest value of a double ~ 2.0E-308
 *      For example:
 *        If 0.0 is not a legitimate number for value, set min = DBL_MIN
 *        If value>=0.0 is legitimate, set min = 0.0
 *        if value<=100., set max = 100.
 *        If no range checking is required, set max = DBL_MAX, min = -DBL_MAX
 *
 *      Certain ascii strings are checked for first (case is insignificant):
 *
 *                 String              Retn_Value ( specified in <limits.h>
 *                 ---------        ---------------------------------------
 *                  FLT_MAX, all          FLT_MAX
 *                  DBL_MAX, max          DBL_MAX
 *                  N/A                  -DBL_MAX
 *                  default               default_value
 *                  small, DBL_MIN        DBL_MIN
 *                  DBL_EPSILON           DBL_EPSILON
 *
 *      A default may be specified. The absence of a
 *      default may also be specified by setting the value of default_value
 *      to NO_DEFAULT_DOUBLE.
 *
 *      If there is an error, *error is set to true. *error isn't touched
 *      if there isn't an error.
 */
{
    double retn_value;
    int    check = false;
    if (defaultVal == NO_DEFAULT_DOUBLE) {
        check = true;
    }
    if (interpret_double(dbl_string, &retn_value, maxVal, minVal, defaultVal)) {
        if (check)
            if (retn_value == NO_DEFAULT_DOUBLE) {
                (void) fprintf(stderr,
                               "ERROR: keyLine_double: Default not allowed\n");
                *error = true;
            }
        if (outofbnds(retn_value, maxVal, minVal)) {
            (void) fprintf(stderr, "ERROR: keyLine_double outofbnds:\n\t\"%s\"\n",
                           dbl_string);
            (void) fprintf(stderr,"\tmax = %e, min = %e\n", maxVal, minVal);
            *error = true;
        }
    } else {
        *error = true;
    }
    return retn_value;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

bool tok_to_boolean(const TOKEN* tokPtr, const int default_value,
                    bool* error)
/*
*      Interprets the first string of a TOKEN structure as a bool.
*      (i.e., true or false).  Returns the interpreted value as the
*       return value.
*      Errors conditions are created if more than one token is found.
*
*      The following character strings are interpreted
*      (case doesn't matter):
*
*      true  = "YES", "true", "T", "Y"
*      false = "NO,   "false", "N", "F"
*      default_value = "DEFAULT" or ""
*
*      A default may be specified on the command line. The absence of a
*      default may also be specified by using the value of NO_DEFAULT_INT.
*      If tokPtr contains no tokens, this routine will try to use the
*      default value.
*
*      If there is an error, *error is set to true. *error isn't touched
*      if there isn't an error.
*/
{
    if (tokPtr->ntokes == 0) {
        return str_to_boolean(DEFAULT_STR, default_value, error);
    } else if (tokPtr->ntokes > 1) {
        (void) fprintf(stderr, "ERROR: tok_to_boolean, ntokes > 1: %s\n",
                       tokPtr->orig_str);
        *error = true;
    }
    return str_to_boolean(tokPtr->tok_ptr[0], default_value, error);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

bool str_to_boolean(const char* string, const int default_value,
                    bool* error)
/*
*      Interprets a stripped character string as a bool value
*      (i.e., true or false). It returns that value as the return value.
*
*      The following character strings are interpreted
*      (case doesn't matter):
*
*      true  = "YES", "true", "T", "Y"
*      false = "NO,   "false", "N", "F"
*      default_value = "DEFAULT"
*
*      A default may be specified on the command line. The absence of a
*      default may also be specified by using the value of NO_DEFAULT_INT.
*
*      If there is an error, *error is set to true. *error isn't touched
*      if there isn't an error.
*/
{
    int retn_value;
    if (interpret_boolean(string, &retn_value, default_value)) {
        if (retn_value == NO_DEFAULT_BOOLEAN) {
            (void) fprintf(stderr,
                           "ERROR: keyLine_boolean: Default not allowed\n");
            *error = true;
        }
    } else {
        *error = true;
    }
    return (retn_value != 0);
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

char* tok_to_string(const TOKEN* tokPtr,  const int maxTok,
                    const int minTok, const char* defaultVal, bool* error)
/*
 *      Interprets the arguments in a TOKEN structure as a string.
 *      It mallocs new space for the string, are returns the pointer to it.
 *      The number of tokens in the string is checked before returning.
 *      The value must be between the maxTok and minTok; it can equal the
 *      max or min value.
 *
 *      Certain ascii strings are checked for first (case is insignificant):
 *
 *                 String              Retn_Value
 *                 ---------        --------------
 *                 default or ""        defaultVal
 *
 *      A default may be specified on the command line. The absence of a
 *      default may also be specified by setting default_value to
 *      NO_DEFAULT_INT.
 *
 *      If there is an error, *error is set to true. *error isn't touched
 *      if there isn't an error.
 */
{
    char* str;
    if (tokPtr->ntokes == 0) {
        str = str_to_string(DEFAULT_STR,      defaultVal, error);
    } else {
        str = str_to_string(tokPtr->orig_str, defaultVal, error);
    }
    if (outofbnds((double) tokPtr->ntokes, (double) maxTok,
                  (double) minTok)) {
        (void) fprintf(stderr, "ERROR: tok_to_String:\n\t\"%s\"\n",
                       tokPtr->orig_str);
        (void) fprintf(stderr,"\tmaxTok = %d, minTok = %d\n", maxTok, minTok);
        *error = true;
    }
    return str;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

char* str_to_string(const char* str, const char* defaultVal,
                    bool* error)
/*
 *      Interprets the argument string as a string.
 *      It mallocs new space for the string, are returns the pointer to it.
 *
 *      Certain ascii strings are checked for first (case is insignificant):
 *
 *                 String              Retn_Value
 *                 ---------        --------------
 *                 default            defaultVal
 *
 *      A default may be specified on the command line. The absence of a
 *      default may also be specified by setting default_value to
 *      NO_DEFAULT_INT.
 *
 *      If there is an error, *error is set to true. *error isn't touched
 *      if there isn't an error.
 */
{
    if (str == NULL) {
        *error = true;
        (void) fprintf(stderr,"ERROR str_to_string: str is uninialized\n");
        return NULL;
    }
    if (strmatch(str, DEFAULT_STR)) {
        if (strmatch(defaultVal, NO_DEFAULT_STR)) {
            *error = true;
            (void) fprintf(stderr,"ERROR str_to_string: no default allowed\n");
            return copy_string(NO_DEFAULT_STR);
        } else {
            return copy_string(defaultVal);
        }
    }
    return copy_string(str);
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

int scan_for_int(FILE* ifp, const char* str, const int maxVal,
                 const int minVal)
/*
 *      Scans the file for a line matching string. Then, interprets
 *      everything after the equals sign as a single integer.
 *      Bounds checking is done on the value before returning.  Value
 *      must be between the max and min; it can equal the max or min value.
 *
 *      Certain ascii strings are checked for first (case is insignificant):
 *
 *                 String              Retn_Value
 *                 ---------        --------------
 *                  INT_MAX, max, all     INT_MAX
 *                  INT_MIN, default      INT_MIN
 *                  N/A, Not Available    INT_MIN
 *
 *      Because this is a fixed format input file routine, errors are
 *      handled by terminally exiting the program.
 */
{
    int retn_value, defaultVal = INT_MIN;
    char* tok_ptr[2];
    char input[MAX_INPUT_STR_LN + 1];
    if (scan_for_line(ifp, str, input, KEY_CHAR, PrintInputFile) < 0) {
        exit(-1);
    }
    if (stokenize(input, DELIMITERS, tok_ptr, 2) == 1) {
        if (interpret_int(tok_ptr[0], &retn_value,
                          maxVal, minVal, defaultVal)) {
            if (outofbnds((double) retn_value, (double) maxVal,
                          (double) minVal)) {
                fprintf(stderr,
                        "ERROR: scan_for_int outofbnds:\n\t\"%s\"\n",str);
                fprintf(stderr,"\tmax = %d, min = %d\n", maxVal, minVal);
                exit(-1);
            } else {
                return retn_value;
            }
        }
    }
    fprintf(stderr, "ERROR while processing line, \"%s\"\n", str);
    exit(-1); /*NOTREACHED*/
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

bool scan_for_boolean(FILE* ifp, const char* string)
/*
 *      Scans the file for a line matching string. Then, interprets
 *      everything after the equals sign as a single boolean value
 *
 *      Because this is a fixed format input file routine, errors are
 *      handled by terminally exiting the program.
 */
{
    int ret_value;
    char* tok_ptr[2];
    char input[MAX_INPUT_STR_LN + 1];
    if (scan_for_line(ifp, string, input, KEY_CHAR, PrintInputFile) < 0) {
        exit(-1);
    }
    if (stokenize(input, DELIMITERS, tok_ptr, 2) == 1) {
        if (interpret_boolean(tok_ptr[0], &ret_value, NO_DEFAULT_BOOLEAN)) {
            if (ret_value == NO_DEFAULT_BOOLEAN) {
                (void) fprintf(stderr, "scan_for_boolean: default not allowed\n");
                exit(-1);
            }
            return (ret_value != 0);
        }
    }
    (void) fprintf(stderr, "scan_for_boolean: ERROR on line \"%s\"\n", string);
    exit(-1);  /*NOTREACHED*/
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double scan_for_double(FILE* ifp, const char* string, const double maxVal,
                       const double minVal)
/*
 *      Scans the file for a line matching string. Then, interprets
 *      everything after the equals sign as a single floating point number.
 *      Bounds checking is then done on the value before returning.  Value
 *      must be between the max and min; it can equal the max or min.
 *
 *      Useful values for bounds:
 *          DBL_MAX = largest legitimate value of a double ~ 2.0E-308
 *         -DBL_MAX = smallest legitimate value of a double ~ 2.0E-308
 *          DBL_EPSILON = smallest value of a double that can be added to one
 *                        and produce a different number. ~ 2.0E-16
 *          DBL_MIN = smallest value of a double ~ 2.0E-308
 *      For example:
 *        If 0.0 is not a legitimate number for value, set min = DBL_MIN
 *        If value>=0.0 is legitimate, set min = 0.0
 *        if value<=100., set max = 100.
 *        If no range checking is required, set max = DBL_MAX, min = -DBL_MAX
 *
 *      Certain ascii strings are checked for first (case is insignificant):
 *
 *                 String              Retn_Value
 *                 ---------        --------------
 *                  FLT_MAX, all          FLT_MAX
 *                  DBL_MAX, max          DBL_MAX
 *                  N/A, default         -DBL_MAX
 *                  small, DBL_MIN        DBL_MIN
 *                  DBL_EPSILON           DBL_EPSILON
 *
 *      Because this is a fixed format input file routine, errors are
 *      handled by terminally exiting the program.
 */
{
    double retn_value;
    char*  tok_ptr[2];
    char input[MAX_INPUT_STR_LN + 1];
    if (scan_for_line(ifp, string, input, KEY_CHAR, true) < 0) {
        exit(-1);
    }
    if (stokenize(input, DELIMITERS, tok_ptr, 2) > 0) {
        if (interpret_double(tok_ptr[0], &retn_value, maxVal, minVal,
                             NO_DEFAULT_DOUBLE)) {
            if (retn_value == NO_DEFAULT_DOUBLE) {
                (void) fprintf(stderr,
                               "ERROR: scan_for_double  has no default\n");
                exit(-1);
            }
            if (outofbnds(retn_value, maxVal, minVal)) {
                (void) fprintf(stderr,
                               "ERROR: scan_for_double outofbnds:\n \"%s = %e\"\n",
                               string, retn_value);
                (void) fprintf(stderr,"\tmax = %e, min = %e\n", maxVal, minVal);
                exit(-1);
            } else {
                return retn_value;
            }
        }
    }
    (void) fprintf(stderr, "ERROR scan_for_double: \"%s\"\n", string);
    exit(-1); /*NOTREACHED*/
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

char* scan_for_string(FILE* ifp, const char* string, const int maxVal,
                      const int minVal)
/*
 *      Scans the file for a line matching string. Then, interprets
 *      everything after the equals sign as string to be returned.
 *      Storage for the resulting string is malloced, and the address
 *      of the string is returned.
 *      The string is returned stripped of leading and trailing white space,
 *      and of comments.
 *
 *      Length checking is then done on the number of characters returned.
 *
 *      Because this is a fixed format input file routine, errors are
 *      handled by terminally exiting the program.
 */
{
    size_t  len;
    char input[MAX_INPUT_STR_LN + 1];
    if (scan_for_line(ifp, string, input, KEY_CHAR, PrintInputFile) < 0) {
        exit(-1);
    }
    len = strlen(input);
    if (outofbnds((double) len, (double) maxVal, (double) minVal)) {
        (void) fprintf(stderr,
                       "ERROR: scan_for_string string length: \"%s = %s\"\n", string, input);
        (void) fprintf(stderr, "\tlength max = %d, min = %d\n", maxVal, minVal);
        exit(-1);
    }
    return copy_string(input);
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

int scan_for_line(FILE* ifp, const char* str, char input[],
                  const char ch_term, const int print_flag)
/*
* Scan the input file (reading in strings according to * 'read_string(ifp,)'
* specifications) until the character pattern in 'string' is matched.
* Returns all of the characters after the termination character in
* a null-character terminated string,
*
* Parameter list:
*
* ifp    == pointer to file "input"
* string == contains string pattern to be matched.
* input  == buffer array to hold characters that are read in.
*                 On output, it contains the return character string
* ch_term== Termination character. When scanning a line of input
* is read until either a newline, the 'ch' termination
* character is read, or the end-of-file is read.
*
* Output:
*     This function returns the number of characters in the string input,
*     excluding the null character.  Error conditions are currently
*     handled by returning with negative return values.
*/
{
    int       retn_value, i;
    bool   found = false;
    char      match_string[MAX_INPUT_STR_LN+1],
              save_input[MAX_INPUT_STR_LN+1];
    static const char* ename = "ERROR scan_for_line: ";

    /*
    *         Error test the input match string
    */

    if (strlen(str) > MAX_INPUT_STR_LN) {
        fprintf(stderr,"%sMatch string is too long:\n\t%s\n",
                ename, str);
        return -1;
    }

    /*
     *        Make it an error to have the comment indicator in a
     *        match string
     */

    for (i = 0; i < (int) strlen(str); i++) {
        if (str[i] == COM_CHAR || str[i] == COM_CHAR2) {
            fprintf(stderr, "%s Comment in match string\n\t%s\n",
                    ename, str);
            return -1;
        }
    }

    /*
    *         Strip the string of leading and trailing white space
    */

    if ((retn_value = strip(strcpy(match_string, str))) <= 0) {
        fprintf(stderr, "%sMatch string is white space: \"%s\"\n",
                ename, str);
        return -1;
    }

    /*
    *          Start the search for the string
    */

    do {
        /*
        *    Read a chunk of text, either up to a newline or to the
        *    character ch_term, from the file pointer, ifp.
        */

        if ((retn_value = read_string(ifp, save_input, ch_term)) < 0) {
            fprintf(stderr,
                    "%sEOF found in input file while searching for:\n", ename);
            fprintf(stderr, "\t\"%s\"\n", match_string);
            return retn_value;
        }

        /*
        *    copy the string just read to the output string, and the
        *    strip it of leading and trailing white space, and comments.
        *    Then, compare the stripped input string with the stripped
        *    match_string
        */

        strcpy(input, save_input);
        if (strip(input) > 0) {
            found =  strmatch(input, match_string);
        }

        /*
        *    If requested print the line, including comments, on standard output.
        *    Use the retn_value from read_string to test whether a \n needs to
        *    be written.
        */

        if (print_flag) {
            if (found) {
                printf("->\t");
            }
            if (retn_value == 0) {
                printf("%s\n", save_input);
            } else {
                printf("%s%c", save_input, ch_term);
            }
        }

        /*
         *   Read and print the rest of the line, if we are in the middle of it.
         */

        if (retn_value > 0) {
            if ((retn_value = read_line(ifp, input, print_flag)) < 0) {
                fprintf(stderr,
                        "ERROR, EOF found in input file while reading line:\n");
                fprintf(stderr, "%s %c\n", str, ch_term);
                return(retn_value);
            }
        } else {
            input[0] = '\0';
        }
    } while (!found);
    return (retn_value);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int read_line(FILE* ifp, char input[], const int print_flag)

/*
*   read_line:
*
*      Reads a line of input.   The line is
*      printed to standard output, if print_flag is true.
*           The line is returned in the character
*      string pointed to by input. Leading and trailing white spaces are
*      stripped from the line.
*          The number of characters, excluding the null character, is
*      returned, except when read_string encounters an error condition
*      (negative return values).  Then, the error condition is returned.
*/
{
    int retn_value;

    /*
    *   read the file up to the next new line, read_string will return
    *   a 0 for a success, and a negative value for failure
    */

    retn_value = read_string(ifp, input, '\n');

    /*
    *   Print out the line before stripping it of comments and white space
    */

    if (print_flag) {
        printf("%s\n", input);
    }

    /*
    *   Strip the return line of comments and leading/trailing white space
    *   Use the function strip to return the number of characters remaining
    */

    if (retn_value == 0) {
        return strip(input);
    }

    /*
    *   If an error condition occurred in read_string, return the error
    *   condition value instead of the character count
    */

    (void) strip(input);
    return retn_value;
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

int read_string(FILE* ifp, char string[], const char ch)
/*
 * This routine reads the standard input until encountering
 * the end-of-file, a newline,  the character 'ch' or until
 * MAX_INPUT_STR_LN characters are read. The inputted characters
 * are read into 'string'.
 *        If an EOF occurs, -1 is returned.
 *        If a line is longer than MAX_INPUT_STR_LN, a -2 is returned
 *        and an error message is written to standard error.
 *        string[] will be returned null-terminated with the
 *        first MAX_INPUT_STR_LN of the line.
 * Upon successful completion with the read terminated by the
 * character 'ch', the number of characters read plus 1 for the
 * null character at the end of the string is returned.  If the
 * read is terminated by '\n', a 0 is returned, even if  ch = '\n'
 *
 *
 * Parameter list:
 *
 * ifp    == pointer to file "input"
 * string == On output, 'string' contains the characters read
 *           from the input stream.  However, the termination character
 *           or the newline character is not included
 * ch     == Additional Termination character. That is, input function
 *           stops when 'ch' or '\n' is read.
 */
{
    int   i = 0, rtn_value, new_ch;

    /*
     *   Read from the file, until termination conditions occur
     *   The order in the while statement is important.
     */
    while ((i < MAX_INPUT_STR_LN)
            && ((new_ch  = getc(ifp)) != ch)
            && (new_ch               != '\n')
            && (new_ch               != EOF)) {
        string[i++] = new_ch;
    }

    /*
     *   Check for termination conditions
     */
    if (new_ch == EOF) {
        rtn_value = -1;
    } else if (i == MAX_INPUT_STR_LN) {
        fprintf(stderr,
                "read_string ERROR: Maxed line character count, %d,"
                " before finding (%c)\n", MAX_INPUT_STR_LN, ch);
        rtn_value = -2;
    } else if (new_ch == '\n') {
        rtn_value =  0;
    } else {
        rtn_value = i + 1;
    }

    /*
     *   Make sure the string is null terminated and return
     */
    string[i] = '\0';
    return rtn_value;
}
/**************************************************************************/

static bool interpret_boolean(const char* token, int* ret_value,
                              const int default_value)
/*
*   This routine interprets a string token to be either true or false
*   and then returns the appropriate answer as an int value in the
*   variable ret_value. It is int because the default value may not
*   be only 0 or 1
*/
{
    /* lower_case (token); */
    if (token[0] == '\0') {
        *ret_value = default_value;
    } else if (strlen(token) == 1) {
        switch (token[0]) {
        case 't':
        case 'y':
            *ret_value = true;
            break;
        case 'f':
        case 'n':
            *ret_value = false;
            break;
        default:
            return false;
        }
    } else {
        if (strmatch(token,"true") || strmatch(token,"yes")) {
            *ret_value = true;
        } else if (strmatch(token,"false") || strmatch(token,"no")) {
            *ret_value = false;
        } else if (strmatch(token,DEFAULT_STR) == 0) {
            *ret_value = default_value;
        } else {
            return false;
        }
    }
    return (true);
}
/**************************************************************************/

static bool interpret_int(const char* token, int* retn_value,
                          const int maxVal,  const int minVal,
                          const int defaultVal)
/*
*   This routine interprets a string token to be an integer
*   and then returns the appropriate answer as an int value in the
*   variable ret_value.
*     Errors are indicated by returning false. Success is indicated
*   by returning true.
*
*      Certain ascii strings are checked for first (case is insignificant):
*
*                 String                 Retn_Value
*                 ---------           --------------
*             INT_MAX, all               INT_MAX
*             INT_MIN                    INT_MIN
*             max                        maxVal
*             min                        minVal
*             default                    defaultVal
*             NULL string                defaultVal
*             N/A, Not_Available         INT_MIN
*/
{
    int      retn;

    /*
     *     Allow a few key ascii phrases in place of an actual int
     */

    /* lower_case(token); */
    if (token[0] == '\0') {
        *retn_value = defaultVal;
    } else if ((strmatch(token,"all")) || strmatch(token,"int_max")) {
        *retn_value = INT_MAX;
    } else if (strmatch(token,"int_min")) {
        *retn_value = INT_MIN;
    } else if (strmatch(token,"max")) {
        *retn_value = maxVal;
    } else if (strmatch(token,"min")) {
        *retn_value = minVal;
    } else if (strmatch(token,DEFAULT_STR)) {
        *retn_value = defaultVal;
    } else if (strmatch(token,"n/a") || strmatch(token,"not_available")) {
        *retn_value = INT_MIN;
    } else  {
        if ((retn = sscanf(token, "%d", retn_value)) != 1) {
            *retn_value = retn;
            return false;
        }
    }
    return true;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

static bool interpret_double(const char* token, double* retn_value,
                             const double maxVal, const double minVal,
                             const double defaultVal)
/*
*   This routine interprets a string token to be a double
*   and then returns the appropriate answer as a double value in the
*   variable retn_value.
*   The function itself returns true if successful or false if unsuccessful
*
*
*      Certain ascii strings are checked for first (case is insignificant):
*
*                 String              Retn_Value
*                 ---------        --------------
*                  FLT_MAX, all          FLT_MAX
*                  FLT_MIN               FLT_MIN
*                  DBL_MAX               DBL_MAX
*                  max                   maxVal
*                  min                   minVal
*                  default               defaultVal
*                  NULL string           defaultVal
*                  N/A                  -DBL_MAX
*                  small, DBL_MIN        DBL_MIN
*                  DBL_EPSILON           DBL_EPSILON
*
*          DBL_MAX = largest legitimate value of a double ~ 2.0E-308
*         -DBL_MAX = smallest legitimate value of a double ~ - 2.0E-308
*          DBL_EPSILON = smallest value of a double that can be added to one
*                        and produce a different number. ~ 2.0E-16
*          DBL_MIN = tiniest value of a double ~ 2.0E-308
*/
{
    int    retn;
    float  retn_float;

    /*
     *     Allow a few key ascii phrases in place of an actual float
     */

    /* lower_case(token); */
    if (token[0] == '\0') {
        *retn_value = defaultVal;
    } else if ((strmatch(token,"all"))      ||  strmatch(token,"flt_max")) {
        *retn_value = FLT_MAX;
    } else if (strmatch(token,"flt_min")) {
        *retn_value = FLT_MIN;
    } else if (strmatch(token,"dbl_max")) {
        *retn_value = DBL_MAX;
    } else if (strmatch(token,"max")) {
        *retn_value = maxVal;
    } else if (strmatch(token,"min")) {
        *retn_value = minVal;
    } else if (strmatch(token,"n/a")) {
        *retn_value = -DBL_MAX;
    } else if (strmatch(token, DEFAULT_STR)) {
        *retn_value = defaultVal;
    } else if (strmatch(token,"small")   || strmatch(token,"dbl_min")) {
        *retn_value = DBL_MIN;
    } else if (strmatch(token,"dbl_epsilon")) {
        *retn_value = DBL_EPSILON;
    } else {
        if ((retn = sscanf(token, "%e", &retn_float)) != 1) {
            *retn_value = (double) retn;
            return false;
        } else {
            *retn_value = (double) retn_float;
        }
    }
    return true;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int strip(char str[])
/*
* This routine strips off blanks and tabs (only leading and trailing
* characters) in 'str'.  On return, it returns the number of
* characters still included in the string (excluding the null character).
*
* Comments are excluded -> All instances of the comment character, '!',
*                          are replaced by '\0' thereby terminating
*                          the string
*
* Parameter list:
*
* str == On output 'str' contains the same characters as on
*        input except the leading and trailing white space and
*        comments have been removed.
*/
{
    int  i = 0, j = 0;
    char ch;

    /*
    *    Quick Returns
    */

    if ((str == NULL) || (str[0] == '\0')) {
        return 0;
    }

    /* Find first non-space character character */

    while (((ch = str[i]) != '\0') && isspace(ch)) {
        i++;
    }

    /*
     * Move real part of str to the front by copying the string
     *   - Comments are handled here, by terminating the copy at the
     *     first comment indicator, and inserting the null character at
     *     that point.
     */

    while ((ch = str[j+i]) != '\0' &&
            (ch != COM_CHAR)           &&
            (ch != COM_CHAR2)) {
        str[j] = ch;
        j++;
    }
    str[j] = '\0';
    j--;

    /* Remove trailing white space by inserting a null character */

    while ((j != -1) && isspace(str[j])) {
        j--;
    }
    j++;
    str[j] = '\0';
    return j;
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/

void lower_case(char str[])
/*
*   lower_case:
*       Translates a string delimited by a NULL character
*       to lower case.  There is no error checking in this version.
*       Relies on stlib function, tolower.
*/
{
    int i;
    for (i = 0; i < (int) strlen(str); i++) {
#   if defined(_INCLUDE_XOPEN_SOURCE) && ! defined(__lint)
        str[i] = _tolower((str[i]));
#   else
        str[i] = tolower(str[i]);
#   endif
    }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

char* TokToStrng(const TOKEN* keyptr)
/*
*    TokToStrng:
*        Mallocs a new character string and copies
*        the tokens character string to it, appending all tokens together
*        into a single string separated by a single space character.
*        It returns the pointer to the new string;
*        The new string should be freed when no longer needed.
*/
{
    int  i;
    if (!keyptr) {
        return NULL;
    }
    if (!keyptr->orig_str) {
        return NULL;
    }
    size_t  iln = strlen(keyptr->orig_str) + 1 + keyptr->ntokes;
    char* fstr = (char*) malloc(iln * sizeof(char));

    char* const* str  = &(keyptr->tok_ptr[0]);
    for (i = 0, fstr[0]= '\0'; i < (keyptr->ntokes - 1); i++, str++) {
        (void) strcat(strcat(fstr, *str), " ");
    }
    return strcat(fstr, *str);
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int stokenize(char* string, const char* delimiters, char* tok_ptr[],
              const int max_tokens)
/*
 * stokenize
 *
 *      This function will break up a string into its respective "tokens".
 *  It returns the number of tokens found. See the strtok(3) man page.
 *
 *  input
 *  ----------
 *    string        - String to be tokenized.  Note, that the string is
 *                    changed by this procedure. Null characters are
 *                    put between each symbol.
 *    delimiters    - String containing a list of delimiters.
 *                    The example below covers 'white space'
 *                      e.g., char *delimiters  = " \t\n";
 *    max_tokens    - Maximum number of tokens to be found
 *
 *  output
 * -----------
 *    tok_ptr    - Vector of pointers to strings, that contain the input
 *                    string's tokens
 */
{
    int i = 0;
    if (string == NULL) {
        tok_ptr[0] = NULL;
        return 0;
    }
    if (strlen(string) == 0) {
        tok_ptr[0] = string;
        return 0;
    }
    if ((tok_ptr[0] = strtok(string, delimiters)) != NULL) {
        do {
            if ((++i) == max_tokens) {
                break;
            }
        } while ((tok_ptr[i] = strtok(NULL, delimiters)) != NULL);
    }
    return i;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static bool outofbnds(const double value, const double maxVal,
                      const double minVal)
/*
 *    This routine checks the bounds of a single double value.
 *    If it is inside or on the bounds, it returns false.
 *    If it is outside the bounds, it returns true.
 */
{
    if ((value <= maxVal) && (value >= minVal)) {
        return false;
    }
    return true;
}
/******************************************************************************
 *
 * strmatch():
 *
 *    This routine checks whether one string is the same as another.
 *    Upper case is transformed into lower case before the comparison is done.
 *    Thus, case doesn't matter in the comparison. However, white space
 *    does matter in this comparison.
 *    If they are,    it returns true
 *    If they aren't, it returns false
 */
bool strmatch(const char* s1, const char* s2)
{
    while (*s1 != '\0') {
#   if defined (_INCLUDE_XOPEN_SOURCE) && ! defined(__lint)
        if (_tolower((*s1++)) != _tolower((*s2++))) {
            return false;
        }
#   else
        if (tolower(*s1) != tolower(*s2)) {
            return false;
        }
        s1++;
        s2++;
#   endif
    }
    if (*s2 != '\0') {
        return false;
    }
    return true;
}

/*****************************************************************************
 *
 * strstrmatch():
 *
 *    This routine checks whether two strings are the same modulo differences
 *    in their white space
 */
bool strstrmatch(const char* s1, const char* s2)
{
    struct TOKEN tmpKeyStruct1, tmpKeyStruct2;
    fillTokStruct(&tmpKeyStruct1, s1);
    fillTokStruct(&tmpKeyStruct2, s2);
    return toktokmatch(&tmpKeyStruct2, &tmpKeyStruct1);
}

/*******************************************************************************
 *
 * strtokmatch():
 *
 *    This routine checks whether a string matches the string contained in
 *    the tokens of a keyLineStr.
 *    White space and case are ignored.
 *    If they are, it returns true
 *    If they aren't, it returns false
 */
bool strtokmatch(const TOKEN* keyptr, const char* s2)
{
    struct TOKEN tmpKeyStruct;
    fillTokStruct(&tmpKeyStruct, s2);
    return toktokmatch(keyptr, &tmpKeyStruct);
}

/**************************************************************************
 *
 * toktokmatch()
 *
 *    This routine checks whether two TOKEN structures contain the
 *    same data up to differences in white space.
 *    Case is ignored as well, as strmatch is called.
 */
bool toktokmatch(const TOKEN* keyptr1, const TOKEN* keyptr2)
{
    int i = keyptr1->ntokes;
    if (i != keyptr2->ntokes) {
        return false;
    }
    while (i > 0) {
        i--;
        if (!strmatch(keyptr1->tok_ptr[i], keyptr2->tok_ptr[i])) {
            return false;
        }
    }
    return true;
}

/**************************************************************************/
/*
 * fillTokStruct()
 *
 *    Fill in a keyLineStruct with a string. Use the defn of white space
 *    at the start of the file to tokenize the string, storing it in the
 *    TOKEN structure.
 */
void fillTokStruct(TOKEN* keyptr1, const char* s2)
{
    if (keyptr1 == NULL) {
        return;
    }
    if (keyptr1->orig_str) {
        free(keyptr1->orig_str);
    }
    if (keyptr1->tok_str) {
        free(keyptr1->tok_str);
    }
    if (s2 == NULL) {
        keyptr1->orig_str = copy_string("");
        keyptr1->tok_str = copy_string("");
        keyptr1->ntokes = 0;
        keyptr1->tok_ptr[0] = keyptr1->orig_str;
        return;
    }
    keyptr1->orig_str = copy_string(s2);
    keyptr1->tok_str  = copy_string(s2);
    keyptr1->ntokes = stokenize(keyptr1->tok_str, DELIMITERS, keyptr1->tok_ptr,
                                MAXTOKENS);
}

/******************************************************************************
 *
 * copyTokStruct():
 *
 *    Copies the information stored in keyptr2 into keyptr1
 */
void copyTokStruct(TOKEN* keyptr1, const TOKEN* keyptr2)
{
    if (keyptr1 == NULL) {
        return;
    }
    if (keyptr2 == NULL) {
        return;
    }
    if (keyptr1->orig_str) {
        free(keyptr1->orig_str);
    }
    if (keyptr1->tok_str) {
        free(keyptr1->tok_str);
    }
    if (keyptr2->orig_str == NULL) {
        keyptr1->orig_str = copy_string("");
        keyptr1->tok_str = copy_string("");
        keyptr1->ntokes = 0;
        keyptr1->tok_ptr[0] = keyptr1->orig_str;
        return;
    }
    keyptr1->orig_str = copy_string(keyptr2->orig_str);
    keyptr1->tok_str  = copy_string(keyptr2->orig_str);
    keyptr1->ntokes = stokenize(keyptr1->tok_str, DELIMITERS, keyptr1->tok_ptr,
                                MAXTOKENS);
}

/******************************************************************************
 *
 * in_char_list():
 *
 *   Finds a match of one string against a list of strings. Returns
 *   the position that the first match occurred.
 *   If no match occurred, returns -1.
 *   The comparisons ignore differences in white space.
 */
int in_char_list(const char* const str1, const char** const list,
                 int num_list)
{
    int i;
    for (i = 0; i < num_list; i++) if (strstrmatch(str1, list[i])) {
            return i;
        }
    return -1;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

char* copy_string(const char* string)
/*
*    copy_string:
*        Mallocs a new character string and copies the old string to it
*
*    NOTE: Memory leak may result if the calling program doesn't free
*          the malloced space
*/
{
    char* new_string;
    new_string = (char*) malloc(strlen(string)+1);
    if (new_string == NULL) {
        (void) fprintf(stderr, "copy_string ERROR: malloc returned NULL");
    } else {
        (void) strcpy(new_string, string);
    }
    return new_string;
}

/******************************************************************************
 *
 * strip_item_from_token ():
 *
 *  Change the token by taking eliminating the iword'th token from the token
 *  structure and reformulating the token expression
 */
void strip_item_from_token(int iword, TOKEN* tok)
{
    if (!tok) {
        return;
    }
    if (iword < 0 || iword > tok->ntokes) {
        return;
    }
#ifdef _MSC_VER
    __w64 size_t ioffset = tok->tok_ptr[iword] - tok->tok_str;
#else
    size_t ioffset = tok->tok_ptr[iword] - tok->tok_str;
#endif
    size_t ilength = strlen(tok->tok_ptr[iword]);
#ifdef _MSC_VER
    __w64 size_t i = ioffset;
    __w64 size_t j = ioffset + ilength;
#else
    size_t i = ioffset;
    size_t j = ioffset + ilength;
#endif
    if (j <= strlen(tok->orig_str)) {
        while (tok->orig_str[j] != '\0') {
            tok->orig_str[i] = tok->orig_str[j];
            i++;
            j++;
        }
        tok->orig_str[i] = '\0';
    }
    strcpy(tok->tok_str, tok->orig_str);
    tok->ntokes = stokenize(tok->tok_str, DELIMITERS, tok->tok_ptr,
                            MAXTOKENS);
}

/*****************************************************************************/
