/*Write a C program to convert specified days into years, weeks and days. Go to the editor
Note: Ignore leap year. 

Test Data :
Number of days : 1329 
Expected Output :
Years: 3 
Weeks: 33 
Days: 3 
*/
#include <stdio.h>
#include <stdlib.h>

int main(){

int num_days = 1329;
int days;
int years;
int weeks;

years = num_days/365;
printf("years: %d\n", years);

weeks = (num_days % 365)/7;
printf("weeks: %d\n", weeks);

days = num_days - ((years*365) + (weeks*7));
printf("days: %d\n", days);
return 0;
}
