#include "bus.hpp"
bus::bus(){
    myCpu.connect(this);
}

bus::~bus(){}

void bus::write(u16 addr, u8 data){
    if(addr >= 0x0000 && addr <= 0xFFFF) ram[addr] = data;
}

u8 bus::read(u16 addr, bool isReadOnly){
    if(addr >= 0x0000 && addr <= 0xFFFF) return ram[addr];
    return 0x00;
}
