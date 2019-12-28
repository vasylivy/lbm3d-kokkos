#include "output.h"

Output::Output() {
}

Output::~Output() {
}

void Output::write_view(const std::string &fn, ScalarField::HostMirror data) {

  char buff[256];

  sprintf(buff, "%s_%03lu.bin", fn.c_str(), frame);

  file = fopen(buff, "wb");

  Double *data_out = data.data();

  if (file != NULL) {
    fwrite(data_out, sizeof(data_out[0]), data.extent(0) * data.extent(1) * data.extent(2), file);
  } else {
    printf("File %s not found\n", fn.c_str());
    exit (EXIT_FAILURE);
  }

  fclose (file);
}

