#include <stdio.h>
#include <string>
#include <map>
#include <inttypes.h>
#define u8 uint8_t
#define u16 uint16_t
class bus;
class cpu{
    public:
        u8  a = 0,
            x = 0,
            y = 0,
            st= 0,
            sp= 0;
        u16 pc= 0x0000;
        enum flags : u8{
            C = 1<<0, // Carry
            Z = 1<<1, // Zero
            I = 1<<2, // Disable Interrupts
            D = 1<<3, // Decimal mode
            B = 1<<4, // Break
            X = 1<<5, // Unused Bit
            V = 1<<6, // Overflow
            N = 1<<7  // Negative
        };

        u8      get(flags x);
        void    assign(bool a, flags x);

        bus *b = nullptr;
        void connect(bus *x){
            b = x;
        };

        u8 read (u16 address);
        void write (u16 address, u8 data);

        struct inst{
            std::string name;
            u8 (cpu::*opcode)(void) = nullptr;
            u8 (cpu::*addrMode)(void) = nullptr;
            u8 cycles = 0;
        };
        inst *lookup;

    private:
        
        void reset();
        void irq();
        void nmi();
        u8 checkCyclesRequired(u8 x);
        void clock();
        bool cmpl();

        u8 opcode = 0x00;
        u8 cycles = 0x00;
        u8 fetch  = 0x00;
        u16 next  = 0x0000;
        u16 nxtr  = 0x0000;
        uint32_t timer = 0;
            
        //Addressing modes:

        u8 imp(); 
        u8 imm();
        u8 zpz();
        u8 zpx();
        u8 zpy();
        u8 rel();
        u8 abs();
        u8 abx();
        u8 aby();
        u8 ind();
        u8 izx();
        u8 izy();

        //Opcodes:

        u8 adc();   u8 nnd();   u8 asl();
        u8 bcc();   u8 bcs();   u8 beq();
        u8 bit();   u8 bmi();   u8 bne();
        u8 bpl();   u8 brk();   u8 bvc();
        u8 bvs();   u8 clc();   u8 cld();
        u8 cli();   u8 clv();   u8 cmp();
        u8 cpx();   u8 cpy();   u8 dec();
        u8 dex();   u8 dey();   u8 eor();
        u8 inc();   u8 inx();   u8 iny();
        u8 jmp();   u8 jsr();   u8 lda();
        u8 ldx();   u8 ldy();   u8 lsr();
        u8 nop();   u8 ora();   u8 pha();
        u8 php();   u8 pla();   u8 plp();
        u8 rol();   u8 ror();   u8 rti();
        u8 rts();   u8 sbc();   u8 sec();
        u8 sed();   u8 sei();   u8 sta();
        u8 stx();   u8 sty();   u8 tax();
        u8 tay();   u8 tsx();   u8 txa();
        u8 txs();   u8 tya();   u8 err();
};