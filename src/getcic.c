#include <stdio.h>
#include <ctype.h>
#include "mytypes.h"
#include "getcic.h"
/****************************************************************
Modify getc() to ignore comments, which begin with # and end with
newline or EOF.  When a comment appears getcic() returns the newline
(or EOF) character only. 
****************************************************************/
int          getcic(FILE * fp)
{
  int          c;

  c = getc(fp);
  if (c == '%')
    do {
      c = getc(fp);
    } while (c != '\n' && c != EOF);
  return (c);
}
/*
 * get a word, delimited by whitespace, from file ifp
 * Ignore comments delimited by %...\n
 */
char        *getwordic(char *buff, int bufsiz, FILE * ifp)
{
  char        *bp;
  int          c;

  bp = buff;
  do {
    c = getcic(ifp);
  } while (isspace(c) && c != EOF);

  if (c == EOF)
    return (NULL);
  while (!isspace(c) && c != EOF) {
    if (--bufsiz < 1)
      break;
    *bp++ = c;
    c = getcic(ifp);
  }
  *bp = '\0';
  return (buff);
}

/* read a floating point number and place its value into x */
int          getrealic(real * x, char *buff, int bufsiz, FILE * ifp)
{
  if (getwordic(buff, bufsiz, ifp) == NULL)
    return (EOF);
  sscanf(buff, RFMT, x);
  return (0);
}

/* read an int and place its value into i */
int          getintic(int *i, char *buff, int bufsiz, FILE * ifp)
{
  if (getwordic(buff, bufsiz, ifp) == NULL)
    return (EOF);
  sscanf(buff, "%d", i);
  return (0);
}
