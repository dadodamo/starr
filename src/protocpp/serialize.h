//
// Created by Daniel Adamovic on 28/10/23.
//

#include<iostream>
#include "../build/proto/ydata.pb.h"
#include "../build/proto/paramdata.pb.h"
#include <fstream>


#ifndef AR_GIBBS_SERIALIZE_H
#define AR_GIBBS_SERIALIZE_H

namespace proto {
    void serialize(sampler_data::samples& sample_stream);
    void serialize_y(y_data::full_y& y_stream);
}

#endif //AR_GIBBS_SERIALIZE_H
