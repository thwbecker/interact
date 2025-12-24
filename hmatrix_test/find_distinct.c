#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/* 

   mainly Claude code generated but tested quite a bit, seems to work

 */
#define TOLERANCE 1e-5
#define INITIAL_CAPACITY 1000
#define MAX_LINE_LENGTH 1024

typedef struct {
  double k;
  double y;
} Entry;

typedef struct {
  Entry *data;
  int size;
  int capacity;
} EntryArray;

void init_array(EntryArray *arr) {
  arr->data = malloc(INITIAL_CAPACITY * sizeof(Entry));
  if (!arr->data) {
    fprintf(stderr, "Memory allocation failed\n");
    exit(1);
  }
  arr->size = 0;
  arr->capacity = INITIAL_CAPACITY;
}

void add_entry(EntryArray *arr, double k, double y) {
  if (arr->size >= arr->capacity) {
    arr->capacity *= 2;
    Entry *new_data = realloc(arr->data, arr->capacity * sizeof(Entry));
    if (!new_data) {
      fprintf(stderr, "Memory reallocation failed\n");
      free(arr->data);
      exit(1);
    }
    arr->data = new_data;
  }
  arr->data[arr->size].k = k;
  arr->data[arr->size].y = y;
  arr->size++;
}

void free_array(EntryArray *arr) {
  free(arr->data);
  arr->data = NULL;
  arr->size = 0;
  arr->capacity = 0;
}

int compare_entries(const void *a, const void *b) {
  Entry *ea = (Entry *)a;
  Entry *eb = (Entry *)b;
  
  if (ea->k < eb->k) return -1;
  if (ea->k > eb->k) return 1;
  return 0;
}

int main(int argc, char *argv[]) {
  int iflag = 0;  // Default value for n filter
  
  if (argc < 2 || argc > 3) {
    fprintf(stderr, "Usage: %s <filename> [iflag]\n", argv[0]);
    fprintf(stderr, "  iflag: integer value to filter on (default: %i)\n",
	    iflag);
    return 1;
  }
  
  // Parse optional iflag argument
  if (argc == 3) {
    iflag = atoi(argv[2]);
  }
  
  FILE *fp = fopen(argv[1], "r");
  if (!fp) {
    fprintf(stderr, "Error opening file: %s\n", argv[1]);
    return 1;
  }
  
  EntryArray entries;
  init_array(&entries);
  
  char line[MAX_LINE_LENGTH];
  int line_num = 0;
  int discarded = 0;
  
  // Read entries line by line
  while (fgets(line, sizeof(line), fp)) {
    line_num++;
    
    double k;
    int n;
    char y_str[64];
    
    // Parse the line - read y as string first to check for NaN
    if (sscanf(line, "%lf %s %d", &k, y_str, &n) == 3) {
      // Check if y_str is "NaN", "nan", or similar
      if (strcasecmp(y_str, "nan") == 0 || 
	  strcasecmp(y_str, "NaN") == 0 ||
	  strcasecmp(y_str, "NAN") == 0) {
	discarded++;
	continue;
      }
      
      // Try to convert y_str to double
      char *endptr;
      double y = strtod(y_str, &endptr);
      
      // Check if conversion was successful and if result is NaN
      if (endptr == y_str || isnan(y)) {
	discarded++;
	continue;
      }
      
      // Only add entries where n == iflag
      if (n == iflag) {
	add_entry(&entries, k, y);
      }
    }
  }
  fclose(fp);
  
  if (discarded > 0) {
    fprintf(stderr, "Discarded %d lines with NaN y values\n", discarded);
  }
  
  if (entries.size == 0) {
    printf("No entries with n == %d found\n", iflag);
    free_array(&entries);
    return 0;
  }
  
  // Sort by k
  qsort(entries.data, entries.size, sizeof(Entry), compare_entries);
  
  // Count distinct y values for each k
  int i = 0;
  while (i < entries.size) {
    double current_k = entries.data[i].k;
    int distinct_y_count = 0;
    
    // Process all entries with this k value
    int j = i;
    while (j < entries.size && entries.data[j].k == current_k) {
      double current_y = entries.data[j].y;
      int is_distinct = 1;
      
      // Check if this y is distinct from previously counted y values
      for (int m = i; m < j; m++) {
	if (fabs(entries.data[m].y - current_y) < TOLERANCE) {
	  is_distinct = 0;
	  break;
	}
      }
      
      if (is_distinct) {
	distinct_y_count++;
      }
      j++;
    }
    
    printf("%10.8f %d\n", current_k, distinct_y_count);
    i = j;
  }
  
  free_array(&entries);
  return 0;
}
