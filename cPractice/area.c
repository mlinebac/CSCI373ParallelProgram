//Write a C program to compute the perimeter and area of a rectangle with a height of 7 inches. and width of 5 inches. Go to the editor
//Expected Output : 
//Perimeter of the rectangle = 24 inches 
//Area of the rectangle = 35 square inches

#include <stdio.h>
#include <stdlib.h>

int main(){

int perimeter;
int area;
int height = 7;
int width = 5;

perimeter = 2*(height + width);
printf("the perimeter is: %d\n", perimeter);
area = height * width;
printf("the area is: %d\n", area);

return 0;
}
