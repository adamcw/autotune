#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

double binomial(int n, int r, double p); 
double odd_errors(int n, double p);
double solve_for_p(int n, double p_change);

int main(int argc, char **argv) {
	float p, pX, pZ, t, time;
	int j, count, cur_d, d, t_check, num_X_changes, num_Z_changes;
	long int num_X_checks, num_Z_checks;

	char *input, *output;
	FILE *fp_in, *fp_out;

	input = NULL;
	output = NULL;

	// Read in command line parameters
	for (j = 1; j < argc; j++) {
		if (strstr(argv[j], "--input")) {
			input = strstr(argv[j], "=");
			if (input) {
				input++;
			} else {
				if (++j < argc) {
					input = argv[j];
				} else {
					printf("--input requires filename. See -h.\n");
					exit(0);
				}
			}
		}
		else if (strstr(argv[j], "--output")) {
			output = strstr(argv[j], "=");
			if (output) {
				output++;
			} else {
				if (++j < argc) {
					output = argv[j];
				} else {
					printf("--output requires filename. See -h.\n");
					exit(0);
				}
			}
		}
		else {
			if (strcmp(argv[j], "-h") && strcmp(argv[j], "--help")) {
				printf("Unknown switch: %s\n", argv[j]);
			}
			printf( "Convert Runs\n\nUsage:\n"
					"  convertruns [--input=INPUT] [--output=OUTPUT]\n\n"
					"Options:\n"
					"  --input=INPUT         A plain-text input file\n"
					"  --output=OUTPUT       A JSON-formatted output file\n");
			exit(0);
		}
	}

	// Open output file or stdout
	fp_out = (output) ? (FILE *)fopen(output, "w") : stdout;
	
	// Open input file or stdin, error if no such file
	if (input) {
		fp_in = (FILE *)fopen(input, "r");
		if (!fp_in) {
			printf("--input: File not found.\n");
			exit(0);
		}
	} else {
		fp_in = stdin;
	}	
	
	cur_d = 0;
	count = 0;
	fprintf(fp_out, "{");
	while (!feof(fp_in)) {
		fscanf(fp_in, "%d %f %d %ld %ld %d %d %f %f\n", &d, &p, &t_check, &num_X_checks, &num_Z_checks, &num_X_changes, &num_Z_changes, &t, &time); 
		
		if (cur_d != d) {
			if (count > 0) {
				fprintf(fp_out, "],\n");
			}
			fprintf(fp_out, "\"%d\": [", d);
			cur_d = d;
		} else if (count > 0) {
			fprintf(fp_out, ",\n");
		}

		if (t_check > 0 && num_X_checks > 0 && num_Z_checks > 0) {
			pX = (float)solve_for_p(t_check, (double)num_X_changes/num_X_checks);
			pZ = (float)solve_for_p(t_check, (double)num_Z_changes/num_Z_checks);
		} else {
			pX = 0;
			pZ = 0;
		}
		fprintf(fp_out, "{\"p\": \"%8.6e\", \"pX\": \"%8.6e\", \"pZ\": \"%8.6e\"}", p, pX, pZ);
		count++;
	}
	fprintf(fp_out, "]}");

	if (fp_out != stdout) {
		fclose(fp_out);
	}
	if (fp_in != stdin) {
		fclose(fp_in);
	}

	return 0;
}

double binomial(int n, int r, double p) {
   int i, j;
   double x, y;

	//printf("binomial n: %d, r: %d, p: %e\n", n, r, p);

   x = 1;
	i = j = 1;
	//printf("x: %e, i: %d, j: %d\n", x, i, j);
	while (i <= r && j <= n-r && x > 1e-18) {
		if (x < 1) {
			y = (n-i+1)*p / i;
			x *= y;
			i++;
			//printf("x: %e, y: %e, i: %d, j: %d\n", x, y, i, j);
		}
		else {
			y = 1-p;
			x *= y;
			j++;
			//printf("x: %e, y: %e, i: %d, j: %d\n", x, y, i, j);
		}
	}

	if (x <= 1e-18) return 0;
	
	while (i <= r && x > 1e-18) {
		if (x < 1) {
			y = (n-i+1)*p / i;
			x *= y;
			i++;
			//printf("x: %e, y: %e, i: %d, j: %d\n", x, y, i, j);
		}
		else {
			y = (n-r+1)*p / r;
			x *= y;
			r--;
			j++;
			//printf("x: %e, y: %e, i: %d, j: %d\n", x, y, i, j);
		}
	}

	if (x <= 1e-18) return 0;
	
	y = pow(1-p, n-r+1-j);
	x *= y;
	//printf("x: %e, y: %e, i: %d, j: %d\n", x, y, i, j);

	//printf("binomial n: %d, r: %d, p: %e, x: %e\n", n, r, p, x);

	return x;
}

double odd_errors(int n, double p) {
   int i;
   double x, y;

	//printf("odd_errors n: %d, p: %e\n", n, p);

   if (n <= 0) return 0;

   x = 0;

	/*
   for (i = 1; i <= n; i+=2) {
      y = binomial(n, i, p);
		x += y;
   }
	*/

	i = n*p;
	if (i%2 == 0) i--;
   for (; i > 0; i-=2) {
      y = binomial(n, i, p);
		if (y <= 1e-18) break;
		x += y;
   }

	i = n*p;
	if (i%2 == 0) i--;
	i += 2;
   for (; i <= n; i+=2) {
      y = binomial(n, i, p);
		if (y <= 1e-18) break;
		x += y;
   }

   return x;
}


double solve_for_p(int n, double p_change) {
   double low;
   double high;
   double epsilon;
   double mid;

	if (p_change == 0) return 0;

	//printf("solve_for_p n: %d, p_change: %e\n", n, p_change);

   epsilon = 1e-18;

   low = epsilon;
   high = 10*p_change/n;
	if (high > 0.5) high = 0.5;

   while (high / low > 1.00001) {
      mid = (low + high) / 2;
      if (odd_errors(n, mid) < p_change) {
         low = mid;
      }
      else {
         high = mid;
      }
   }

   return (low + high) / 2;
}

