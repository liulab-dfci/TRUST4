#ifndef _MOURISL_BARCODETRANSLATOR
#define _MOURISL_BARCODETRANSLATOR

#include <cstring>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include <zlib.h>

#include "khash.h"

KHASH_MAP_INIT_STR(barcodeHash, char *)

// The class for handling barcode convertion.
class BarcodeTranslator {
 public:
  BarcodeTranslator() {
    barcode_translate_table_ = NULL;
    from_bc_length_ = -1;
  }

  ~BarcodeTranslator() {
    if (barcode_translate_table_ != NULL) {
      khiter_t k;
      for (k = kh_begin(barcode_translate_table_);
           k != kh_end(barcode_translate_table_); ++k) {
        if (kh_exist(barcode_translate_table_, k))
          free(kh_value(barcode_translate_table_, k));
      }
      kh_destroy(barcodeHash, barcode_translate_table_);

			int size = key_pointers_.size() ;
			for (int i = 0 ; i < size ; ++i)
				free(key_pointers_[i]) ;
    }
  }

  void SetTranslateTable(const std::string &file) {
    barcode_translate_table_ = kh_init(barcodeHash);

    gzFile barcode_translate_file = gzopen(file.c_str(), "r");
    const uint32_t line_buffer_size = 512;
    char file_line[line_buffer_size];
    while (gzgets(barcode_translate_file, file_line, line_buffer_size) != NULL) {
      int line_len = strlen(file_line);
      if (file_line[line_len - 1] == '\n') {
        file_line[line_len - 1] = '\0';
      }
      std::string tmp_string(file_line);
      ProcessTranslateFileLine(tmp_string);
    }
    gzclose(barcode_translate_file) ;
  }

  std::string Translate(char *bc, uint32_t bc_length) {
    if (barcode_translate_table_ == NULL) {
			std::string tmp(bc) ;
      return tmp ;
    }

    std::string ret;
    uint64_t i;
    for (i = 0; i < bc_length / from_bc_length_; ++i) {
			std::string bc_from(bc + i * from_bc_length_, from_bc_length_) ;
			khiter_t barcode_translate_table_iter =
          kh_get(barcodeHash, barcode_translate_table_, bc_from.c_str());
      if (barcode_translate_table_iter == kh_end(barcode_translate_table_)) {
        std::cerr << "Barcode " << bc_from << " does not exist in the translation table."
                  << std::endl;
        exit(-1);
      }
      std::string bc_to(
          kh_value(barcode_translate_table_, barcode_translate_table_iter));
      if (i == 0) {
        ret = bc_to;
      } else {
        ret += "-" + bc_to;
      }
    }
    return ret;
  }

	bool isSet()
	{
		return barcode_translate_table_ != NULL ;
	}
 private:
  khash_t(barcodeHash) * barcode_translate_table_;
  int from_bc_length_;
	std::vector<char *> key_pointers_ ;

  void ProcessTranslateFileLine(std::string &line) {
    int i;
    int len = line.length();
    std::string from, to;
    for (i = 0; i < len; ++i) {
      if (line[i] == ',' || line[i] == '\t' || line[i] == ' ') break;
    }

    to = line.substr(0, i);
    from = line.substr(i + 1, len - i - 1);
    from_bc_length_ = len - i - 1;
    int khash_return_code;
		char *tmpStr = strdup(from.c_str()) ;
		key_pointers_.push_back(tmpStr) ;
    khiter_t barcode_translate_table_iter = kh_put(
        barcodeHash, barcode_translate_table_, tmpStr, &khash_return_code);
    kh_value(barcode_translate_table_, barcode_translate_table_iter) =
        strdup(to.c_str());
  }
};

#endif
