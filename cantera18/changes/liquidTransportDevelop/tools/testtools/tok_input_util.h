/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: tok_input_util.h,v $
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 * Symbolic $Name:  $
 *
 * $Id$
 * $Source: /cvsroot/cantera/cantera/tools/testtools/tok_input_util.h,v $
 *
 *
 *====================================================================*/

#ifndef CKTOOLS_TOK_INPUT_UTIL_H
#define CKTOOLS_TOK_INPUT_UTIL_H

#include "stdio.h"

#ifndef MAX_INPUT_STR_LN
#define MAX_INPUT_STR_LN 16100
#endif
#ifndef MAX_TOKEN_STR_LN
#define MAX_TOKEN_STR_LN 255
#endif
#ifndef MAXTOKENS
#define MAXTOKENS   20
#endif

#ifndef TRUE
#  define TRUE  1
#  define FALSE 0
#endif

#ifndef BOOLEAN
# define BOOLEAN int
#endif

/**************************************************************************/
/* TOKEN structure:  
 *  This structure is used to parse strings. The original string is
 *  tokenized into a set of tokens via separation wrt white space. 
 *  Both the tokens and the original 
 *  string are storred within the structure
 */

struct TOKEN {
  char *orig_str;
  char *tok_str;
  char *tok_ptr[MAXTOKENS];
  int   ntokes;
  TOKEN(void);
  ~TOKEN();
  TOKEN(const char *str);
};


#define NO_DEFAULT_INT     -68361
#define NO_DEFAULT_BOOLEAN NO_DEFAULT_INT
#define NO_DEFAULT_DOUBLE  -1.241E+11
#define NO_DEFAULT_STR    "NO_DEFAULT"


/**************************************************************************/
/*           Prototypes for routines tok_input_util.c                     */
/**************************************************************************/

extern BOOLEAN get_next_keyLine(FILE *, TOKEN *, TOKEN *);
extern int     tok_to_int    (const TOKEN *, const int, const int, 
                              const int,    BOOLEAN *);
extern int     str_to_int    (const char *,    const int,    const int, 
                              const int,    BOOLEAN *);
extern double  tok_to_double (const TOKEN *, const double, const double, 
                              const double, BOOLEAN *);
extern double  str_to_double (const char *,    const double, const double, 
                              const double, BOOLEAN *);
extern BOOLEAN tok_to_boolean(const TOKEN *, const int, BOOLEAN *);
extern BOOLEAN str_to_boolean(const char *, const int, BOOLEAN *);
extern char   *tok_to_string (const TOKEN *,  const int,  const int,
                              const char *, BOOLEAN *);
extern char   *str_to_string (const char *, const char *, BOOLEAN *);
extern int    scan_for_int(FILE *, const char *, const int, const int);
extern double scan_for_double(FILE *,  const char *, const double, 
                              const double);
extern char  *scan_for_string(FILE *,  const char *, const int, const int);
extern BOOLEAN scan_for_boolean(FILE *, const char *);
extern int    scan_for_line(FILE *, const char *, char [], const char, 
                            const int);
extern int  read_line(FILE *, char [], const int);
extern int  read_string(FILE *, char [], const char);
extern int strip(char []);
extern void lower_case(char []);
extern char *TokToStrng(const TOKEN *);
extern int  stokenize(char *, const char *, char *[], const int);
extern BOOLEAN strmatch(const char *, const char*);
extern BOOLEAN strstrmatch(const char *, const char*);
extern BOOLEAN strtokmatch(const TOKEN *, const char *);
extern BOOLEAN toktokmatch(const TOKEN *, const TOKEN *);
extern void fillTokStruct(TOKEN *, const char *);
extern void copyTokStruct(TOKEN *, const TOKEN *);
extern int in_char_list(const char * const, const char ** const, int);
extern char *copy_string(const char *);
extern void strip_item_from_token(int, TOKEN *);
/**************************************************************************/
#endif /* END OF TOK_INPUT_UTIL_H */
/**************************************************************************/

