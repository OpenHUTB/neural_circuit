/* Cradd - add carriage return before all linefeeds. */

#include <stdio.h>
#include <stdlib.h>

int main(void)
{
 int c;

 while ((c = getchar()) != EOF)
  {
   if (c == '\n') putchar ('\r');
   if (c >= ' ' || c == 0177 || c == '\t' || c == '\n')
    putchar(c);
  }
 exit (0);
}

