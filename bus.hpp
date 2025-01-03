#include <inttypes.h>
#include "cpu.hpp"
#define u8 uint8_t
#define u16 uint16_t
class bus{
    public:
        bus(){};
        ~bus(){};

        cpu myCpu;
        u8 ram[0x10000];

        void write(u16 address, u8 data);
        u8 read(u16 address, bool isReadOnly = false);
};

// 0011 1010 1010 1001 1111 1010 1010 1010