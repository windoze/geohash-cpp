//
//  geohash.cpp
//
//  Created by Xu, Chen on 14/10/31.
//  Copyright (c) 2014 0d0a.com. All rights reserved.
//

#include <cmath>
#include <algorithm>
#include "geohash.hpp"

////////////////////////////////////////////////////////////////////////////////
// geolocation
////////////////////////////////////////////////////////////////////////////////

constexpr long double EARTH_RADIUS=6371.393;
constexpr long double LONGITUDE_CIRCLE=EARTH_RADIUS*2*M_PI;

inline long double radians(long double d) {
    return d * M_PI / 180;
}

inline long double degrees(long double d) {
    return d * 180 / M_PI;
}

inline long double to180(long double d) {
    return d-(std::llround(d)/360*360);
}

inline long double latitude_circle(long double lat) {
    return EARTH_RADIUS*std::cos(radians(lat)) * 2 * M_PI;
}

inline long double longitude_span(long double lat, long double lon1, long double lon2) {
    return radians(std::abs(to180(lon1-lon2)))*latitude_circle(lat)/2/M_PI;
}

inline long double latitide_span(long double lat1, long double lat2) {
    return radians(std::abs(to180(lat1-lat2)))*LONGITUDE_CIRCLE/2/M_PI;
}

double distance(const geolocation &l, const geolocation &r) {
    long double lat1 = radians(l.latitude);
    long double lat2 = radians(r.latitude);
    long double dlat = radians(r.latitude - l.latitude);
    long double dlon = radians(r.longitude - l.longitude);
    long double a = sin(dlat/2.0) * sin(dlat/2.0) + cos(lat1) * cos(lat2) * sin(dlon/2.0) * sin(dlon/2.0);
    long double c = 2.0 * atan2(sqrt(a), sqrt(1.0-a));
    
    return EARTH_RADIUS * c;
}

////////////////////////////////////////////////////////////////////////////////
// bounding_box
////////////////////////////////////////////////////////////////////////////////

bounding_box::bounding_box(geolocation l, double distance) {
    long double latitude_range=degrees(distance/LONGITUDE_CIRCLE);
    long double max_lat=std::max(std::abs(l.latitude-latitude_range), std::abs(l.latitude+latitude_range));
    long double longitude_range=degrees(distance/latitude_circle(max_lat));
    *this=bounding_box(l.latitude-latitude_range,
                       l.latitude+latitude_range,
                       l.longitude-longitude_range,
                       l.longitude+longitude_range);
}

double bounding_box::min_span() const {
    // Returns the shortest side of the box
    return std::min(latitide_span(min_lat, max_lat),
                    std::min(longitude_span(min_lat, min_lon, max_lon),
                             longitude_span(max_lat, min_lon, max_lon)));
}

////////////////////////////////////////////////////////////////////////////////
// geohash
////////////////////////////////////////////////////////////////////////////////

static const char base32_codes[] = {
    '0', '1', '2', '3', '4', '5', '6', '7',
    '8', '9', 'b', 'c', 'd', 'e', 'f', 'g',
    'h', 'j', 'k', 'm', 'n', 'p', 'q', 'r',
    's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
};

static const int base32_indexes[]={
     0,  1,  2,  3,  4,  5,  6,  7, // 30-37, '0'..'7'
     8,  9, -1, -1, -1, -1, -1, -1, // 38-2F, '8','9'
    -1, -1, 10, 11, 12, 13, 14, 15, // 40-47, 'B'..'G'
    16, -1, 17, 18, -1, 19, 20, -1, // 48-4F, 'H','J','K','M','N'
    21, 22, 23, 24, 25, 26, 27, 28, // 50-57, 'P'..'W'
    29, 30, 31, -1, -1, -1, -1, -1, // 58-5F, 'X'..'Z'
    -1, -1, 10, 11, 12, 13, 14, 15, // 60-67, 'b'..'g'
    16, -1, 17, 18, -1, 19, 20, -1, // 68-6F, 'h','j','k','m','n'
    21, 22, 23, 24, 25, 26, 27, 28, // 70-77, 'p'..'w'
    29, 30, 31,                     // 78-7A, 'x'..'z'
};

inline long double calc_height_degrees(size_t n) {
    return 180/std::pow((long double)(2.0), 2.5*n+((n%2==0)?0:-0.5));
}

inline long double calc_width_degrees(size_t n) {
    return 180/std::pow((long double)(2.0), 2.5*n+((n%2==0)?-1:-0.5));
};

static const long double height_degrees[]={
    calc_height_degrees(0),
    calc_height_degrees(1),
    calc_height_degrees(2),
    calc_height_degrees(3),
    calc_height_degrees(4),
    calc_height_degrees(5),
    calc_height_degrees(6),
    calc_height_degrees(7),
    calc_height_degrees(8),
    calc_height_degrees(9),
    calc_height_degrees(10),
    calc_height_degrees(11),
    calc_height_degrees(12),
};

static const long double width_degrees[]={
    calc_width_degrees(0),
    calc_width_degrees(1),
    calc_width_degrees(2),
    calc_width_degrees(3),
    calc_width_degrees(4),
    calc_width_degrees(5),
    calc_width_degrees(6),
    calc_width_degrees(7),
    calc_width_degrees(8),
    calc_width_degrees(9),
    calc_width_degrees(10),
    calc_width_degrees(11),
    calc_width_degrees(12),
};

binary_hash::binary_hash(const std::string &bitstring)
: bits(0)
, precision(bitstring.size())
{
    for(auto c : bitstring) {
        if (c=='1') {
            bits = (bits<<1) | 1;
        } else {
            bits = bits<<1;
        }
    }
}

binary_hash::binary_hash(geolocation l, double dist) {
    bounding_box box(l, dist);
    for (size_t i = MAX_BINHASH_LENGTH; i >= 1; i--) {
        binary_hash hash=binary_encode(box.bottom_left(), i);
        if (decode(hash).contains(box.top_right())) {
            *this=hash;
            return;
        }
    }
}

std::string binary_hash::to_string() const {
    std::string output(' ', precision);
    uint64_t b=bits;
    
    for (int i=0; i<precision; i++) {
        output[i]=(b & (uint64_t(-1))) ? '1' : '0';
        b<<=1;
    }
    
    return output;
}

binary_hash binary_encode(geolocation l, size_t precision) {
    // bbox for the lat/lon + errors/ranges
    bounding_box bbox{ -90, 90, -180, 180 };
    bool islon = true;
    
    binary_hash output;
    
    while(output.size() < precision) {
        if (islon) {
            if(l.longitude > bbox.lon_center()) {
                output.push_back(true);
                bbox.min_lon=bbox.lon_center();
            } else {
                output.push_back(false);
                bbox.max_lon=bbox.lon_center();
            }
        } else {
            if(l.latitude > bbox.lat_center() ) {
                output.push_back(true);
                bbox.min_lat = bbox.lat_center();
            } else {
                output.push_back(false);
                bbox.max_lat = bbox.lat_center();
            }
        }
        islon = !islon;
    }
    return output;
}

std::string encode(geolocation l, size_t precision) {
    // DecodedBBox for the lat/lon + errors
    bounding_box bbox{ -90, 90, -180, 180 };
    bool islon = true;
    int num_bits = 0;
    int hash_index = 0;
    
    // Pre-Allocate the hash string
    std::string output(precision, ' ');
    size_t output_length = 0;
    
    while(output_length < precision) {
        if (islon) {
            if(l.longitude > bbox.lon_center()) {
                hash_index = (hash_index << 1) + 1;
                bbox.min_lon=bbox.lon_center();
            } else {
                hash_index = (hash_index << 1) + 0;
                bbox.max_lon=bbox.lon_center();
            }
        } else {
            if(l.latitude > bbox.lat_center() ) {
                hash_index = (hash_index << 1) + 1;
                bbox.min_lat = bbox.lat_center();
            } else {
                hash_index = (hash_index << 1) + 0;
                bbox.max_lat = bbox.lat_center();
            }
        }
        islon = !islon;
        
        ++num_bits;
        if (5 == num_bits) {
            output[output_length] = base32_codes[hash_index];
            output_length++;
            num_bits = 0;
            hash_index = 0;
        }
    }
    return output;
}

bounding_box decode(const binary_hash &hash) {
    // bbox for the lat/lon + errors/ranges
    bounding_box output{ -90, 90, -180, 180 };
    
    bool islon = true;
    
    for(size_t i=0; i<hash.size(); i++) {
        bool bit = hash.test(i);
        if (islon) {
            if(bit) {
                output.min_lon = output.lon_center();
            } else {
                output.max_lon = output.lon_center();
            }
        } else {
            if(bit) {
                output.min_lon = output.lat_center();
            } else {
                output.max_lon = output.lat_center();
            }
        }
        islon = !islon;
    }
    return output;
}

bounding_box decode(const std::string &hash) {
    bounding_box output{ -90, 90, -180, 180 };
    
    bool islon = true;
    
    for(auto &c : hash) {
        int char_index = base32_indexes[c-48];
        if (char_index<'0' || char_index>'z') {
            throw std::invalid_argument("Invalid geohash");
        }
        
        for (int bits = 4; bits >= 0; --bits) {
            int bit = (char_index >> bits) & 1;
            if (islon) {
                if(bit == 1) {
                    output.min_lon = output.lon_center();
                } else {
                    output.max_lon = output.lon_center();
                }
            } else {
                if(bit == 1) {
                    output.min_lon = output.lat_center();
                } else {
                    output.max_lon = output.lat_center();
                }
            }
            islon = !islon;
        }
    }
    return output;
}

bool hash_contains(const std::string &hash, geolocation l) {
    return decode(hash).contains(l);
}

size_t hash_precision(geolocation l, double dist) {
    for (size_t i = MAX_GEOHASH_LENGTH; i >= 1; i--) {
        bounding_box box=decode(encode(l, i));
        if (box.min_span()>=dist*2) {
            return i;
        }
    }
    return 0;
}

std::string base_hash(geolocation l, double dist) {
    for (size_t i = MAX_GEOHASH_LENGTH; i >= 1; i--) {
        std::string hash=encode(l, i);
        if (decode(hash).min_span()>=dist*2) {
            return hash;
        }
    }
    return "";
}

size_t binary_hash_precision(geolocation l, double dist) {
    for (size_t i = MAX_BINHASH_LENGTH; i >= 1; i--) {
        bounding_box box=decode(binary_encode(l, i));
        if (box.min_span()>=dist*2) {
            return i;
        }
    }
    return 0;
}
