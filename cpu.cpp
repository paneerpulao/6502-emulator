#include "bus.hpp"
#include "cpu.hpp"

cpu::cpu(){
    using a = cpu;
    lookup = 
    {
        { "brk", &a::brk, &a::imm, 7 },{ "ora", &a::ora, &a::izx, 6 },{ "???", &a::err, &a::imp, 2 },{ "???", &a::err, &a::imp, 8 },{ "???", &a::nop, &a::imp, 3 },{ "ora", &a::ora, &a::zpz, 3 },{ "asl", &a::asl, &a::zpz, 5 },{ "???", &a::err, &a::imp, 5 },{ "php", &a::php, &a::imp, 3 },{ "ora", &a::ora, &a::imm, 2 },{ "asl", &a::asl, &a::imp, 2 },{ "???", &a::err, &a::imp, 2 },{ "???", &a::nop, &a::imp, 4 },{ "ora", &a::ora, &a::abs, 4 },{ "asl", &a::asl, &a::abs, 6 },{ "???", &a::err, &a::imp, 6 },
        { "bpl", &a::bpl, &a::rel, 2 },{ "ora", &a::ora, &a::izy, 5 },{ "???", &a::err, &a::imp, 2 },{ "???", &a::err, &a::imp, 8 },{ "???", &a::nop, &a::imp, 4 },{ "ora", &a::ora, &a::zpx, 4 },{ "asl", &a::asl, &a::zpx, 6 },{ "???", &a::err, &a::imp, 6 },{ "clc", &a::clc, &a::imp, 2 },{ "ora", &a::ora, &a::aby, 4 },{ "???", &a::nop, &a::imp, 2 },{ "???", &a::err, &a::imp, 7 },{ "???", &a::nop, &a::imp, 4 },{ "ora", &a::ora, &a::abx, 4 },{ "asl", &a::asl, &a::abx, 7 },{ "???", &a::err, &a::imp, 7 },
        { "jsr", &a::jsr, &a::abs, 6 },{ "nnd", &a::nnd, &a::izx, 6 },{ "???", &a::err, &a::imp, 2 },{ "???", &a::err, &a::imp, 8 },{ "bit", &a::bit, &a::zpz, 3 },{ "nnd", &a::nnd, &a::zpz, 3 },{ "rol", &a::rol, &a::zpz, 5 },{ "???", &a::err, &a::imp, 5 },{ "plp", &a::plp, &a::imp, 4 },{ "nnd", &a::nnd, &a::imm, 2 },{ "rol", &a::rol, &a::imp, 2 },{ "???", &a::err, &a::imp, 2 },{ "bit", &a::bit, &a::abs, 4 },{ "nnd", &a::nnd, &a::abs, 4 },{ "rol", &a::rol, &a::abs, 6 },{ "???", &a::err, &a::imp, 6 },
        { "bmi", &a::bmi, &a::rel, 2 },{ "nnd", &a::nnd, &a::izy, 5 },{ "???", &a::err, &a::imp, 2 },{ "???", &a::err, &a::imp, 8 },{ "???", &a::nop, &a::imp, 4 },{ "nnd", &a::nnd, &a::zpx, 4 },{ "rol", &a::rol, &a::zpx, 6 },{ "???", &a::err, &a::imp, 6 },{ "sec", &a::sec, &a::imp, 2 },{ "nnd", &a::nnd, &a::aby, 4 },{ "???", &a::nop, &a::imp, 2 },{ "???", &a::err, &a::imp, 7 },{ "???", &a::nop, &a::imp, 4 },{ "nnd", &a::nnd, &a::abx, 4 },{ "rol", &a::rol, &a::abx, 7 },{ "???", &a::err, &a::imp, 7 },
        { "rti", &a::rti, &a::imp, 6 },{ "eor", &a::eor, &a::izx, 6 },{ "???", &a::err, &a::imp, 2 },{ "???", &a::err, &a::imp, 8 },{ "???", &a::nop, &a::imp, 3 },{ "eor", &a::eor, &a::zpz, 3 },{ "lsr", &a::lsr, &a::zpz, 5 },{ "???", &a::err, &a::imp, 5 },{ "pha", &a::pha, &a::imp, 3 },{ "eor", &a::eor, &a::imm, 2 },{ "lsr", &a::lsr, &a::imp, 2 },{ "???", &a::err, &a::imp, 2 },{ "jmp", &a::jmp, &a::abs, 3 },{ "eor", &a::eor, &a::abs, 4 },{ "lsr", &a::lsr, &a::abs, 6 },{ "???", &a::err, &a::imp, 6 },
        { "bvc", &a::bvc, &a::rel, 2 },{ "eor", &a::eor, &a::izy, 5 },{ "???", &a::err, &a::imp, 2 },{ "???", &a::err, &a::imp, 8 },{ "???", &a::nop, &a::imp, 4 },{ "eor", &a::eor, &a::zpx, 4 },{ "lsr", &a::lsr, &a::zpx, 6 },{ "???", &a::err, &a::imp, 6 },{ "cli", &a::cli, &a::imp, 2 },{ "eor", &a::eor, &a::aby, 4 },{ "???", &a::nop, &a::imp, 2 },{ "???", &a::err, &a::imp, 7 },{ "???", &a::nop, &a::imp, 4 },{ "eor", &a::eor, &a::abx, 4 },{ "lsr", &a::lsr, &a::abx, 7 },{ "???", &a::err, &a::imp, 7 },
        { "rts", &a::rts, &a::imp, 6 },{ "adc", &a::adc, &a::izx, 6 },{ "???", &a::err, &a::imp, 2 },{ "???", &a::err, &a::imp, 8 },{ "???", &a::nop, &a::imp, 3 },{ "adc", &a::adc, &a::zpz, 3 },{ "ror", &a::ror, &a::zpz, 5 },{ "???", &a::err, &a::imp, 5 },{ "pla", &a::pla, &a::imp, 4 },{ "adc", &a::adc, &a::imm, 2 },{ "ror", &a::ror, &a::imp, 2 },{ "???", &a::err, &a::imp, 2 },{ "jmp", &a::jmp, &a::ind, 5 },{ "adc", &a::adc, &a::abs, 4 },{ "ror", &a::ror, &a::abs, 6 },{ "???", &a::err, &a::imp, 6 },
        { "bvs", &a::bvs, &a::rel, 2 },{ "adc", &a::adc, &a::izy, 5 },{ "???", &a::err, &a::imp, 2 },{ "???", &a::err, &a::imp, 8 },{ "???", &a::nop, &a::imp, 4 },{ "adc", &a::adc, &a::zpx, 4 },{ "ror", &a::ror, &a::zpx, 6 },{ "???", &a::err, &a::imp, 6 },{ "sei", &a::sei, &a::imp, 2 },{ "adc", &a::adc, &a::aby, 4 },{ "???", &a::nop, &a::imp, 2 },{ "???", &a::err, &a::imp, 7 },{ "???", &a::nop, &a::imp, 4 },{ "adc", &a::adc, &a::abx, 4 },{ "ror", &a::ror, &a::abx, 7 },{ "???", &a::err, &a::imp, 7 },
        { "???", &a::nop, &a::imp, 2 },{ "sta", &a::sta, &a::izx, 6 },{ "???", &a::nop, &a::imp, 2 },{ "???", &a::err, &a::imp, 6 },{ "sty", &a::sty, &a::zpz, 3 },{ "sta", &a::sta, &a::zpz, 3 },{ "stx", &a::stx, &a::zpz, 3 },{ "???", &a::err, &a::imp, 3 },{ "dey", &a::dey, &a::imp, 2 },{ "???", &a::nop, &a::imp, 2 },{ "txa", &a::txa, &a::imp, 2 },{ "???", &a::err, &a::imp, 2 },{ "sty", &a::sty, &a::abs, 4 },{ "sta", &a::sta, &a::abs, 4 },{ "stx", &a::stx, &a::abs, 4 },{ "???", &a::err, &a::imp, 4 },
        { "bcc", &a::bcc, &a::rel, 2 },{ "sta", &a::sta, &a::izy, 6 },{ "???", &a::err, &a::imp, 2 },{ "???", &a::err, &a::imp, 6 },{ "sty", &a::sty, &a::zpx, 4 },{ "sta", &a::sta, &a::zpx, 4 },{ "stx", &a::stx, &a::zpy, 4 },{ "???", &a::err, &a::imp, 4 },{ "tya", &a::tya, &a::imp, 2 },{ "sta", &a::sta, &a::aby, 5 },{ "txs", &a::txs, &a::imp, 2 },{ "???", &a::err, &a::imp, 5 },{ "???", &a::nop, &a::imp, 5 },{ "sta", &a::sta, &a::abx, 5 },{ "???", &a::err, &a::imp, 5 },{ "???", &a::err, &a::imp, 5 },
        { "ldy", &a::ldy, &a::imm, 2 },{ "lda", &a::lda, &a::izx, 6 },{ "ldx", &a::ldx, &a::imm, 2 },{ "???", &a::err, &a::imp, 6 },{ "ldy", &a::ldy, &a::zpz, 3 },{ "lda", &a::lda, &a::zpz, 3 },{ "ldx", &a::ldx, &a::zpz, 3 },{ "???", &a::err, &a::imp, 3 },{ "tay", &a::tay, &a::imp, 2 },{ "lda", &a::lda, &a::imm, 2 },{ "tax", &a::tax, &a::imp, 2 },{ "???", &a::err, &a::imp, 2 },{ "ldy", &a::ldy, &a::abs, 4 },{ "lda", &a::lda, &a::abs, 4 },{ "ldx", &a::ldx, &a::abs, 4 },{ "???", &a::err, &a::imp, 4 },
        { "bcs", &a::bcs, &a::rel, 2 },{ "lda", &a::lda, &a::izy, 5 },{ "???", &a::err, &a::imp, 2 },{ "???", &a::err, &a::imp, 5 },{ "ldy", &a::ldy, &a::zpx, 4 },{ "lda", &a::lda, &a::zpx, 4 },{ "ldx", &a::ldx, &a::zpy, 4 },{ "???", &a::err, &a::imp, 4 },{ "clv", &a::clv, &a::imp, 2 },{ "lda", &a::lda, &a::aby, 4 },{ "tsx", &a::tsx, &a::imp, 2 },{ "???", &a::err, &a::imp, 4 },{ "ldy", &a::ldy, &a::abx, 4 },{ "lda", &a::lda, &a::abx, 4 },{ "ldx", &a::ldx, &a::aby, 4 },{ "???", &a::err, &a::imp, 4 },
        { "cpy", &a::cpy, &a::imm, 2 },{ "cmp", &a::cmp, &a::izx, 6 },{ "???", &a::nop, &a::imp, 2 },{ "???", &a::err, &a::imp, 8 },{ "cpy", &a::cpy, &a::zpz, 3 },{ "cmp", &a::cmp, &a::zpz, 3 },{ "dec", &a::dec, &a::zpz, 5 },{ "???", &a::err, &a::imp, 5 },{ "iny", &a::iny, &a::imp, 2 },{ "cmp", &a::cmp, &a::imm, 2 },{ "dex", &a::dex, &a::imp, 2 },{ "???", &a::err, &a::imp, 2 },{ "cpy", &a::cpy, &a::abs, 4 },{ "cmp", &a::cmp, &a::abs, 4 },{ "dec", &a::dec, &a::abs, 6 },{ "???", &a::err, &a::imp, 6 },
        { "bne", &a::bne, &a::rel, 2 },{ "cmp", &a::cmp, &a::izy, 5 },{ "???", &a::err, &a::imp, 2 },{ "???", &a::err, &a::imp, 8 },{ "???", &a::nop, &a::imp, 4 },{ "cmp", &a::cmp, &a::zpx, 4 },{ "dec", &a::dec, &a::zpx, 6 },{ "???", &a::err, &a::imp, 6 },{ "cld", &a::cld, &a::imp, 2 },{ "cmp", &a::cmp, &a::aby, 4 },{ "nop", &a::nop, &a::imp, 2 },{ "???", &a::err, &a::imp, 7 },{ "???", &a::nop, &a::imp, 4 },{ "cmp", &a::cmp, &a::abx, 4 },{ "dec", &a::dec, &a::abx, 7 },{ "???", &a::err, &a::imp, 7 },
        { "cpx", &a::cpx, &a::imm, 2 },{ "sbc", &a::sbc, &a::izx, 6 },{ "???", &a::nop, &a::imp, 2 },{ "???", &a::err, &a::imp, 8 },{ "cpx", &a::cpx, &a::zpz, 3 },{ "sbc", &a::sbc, &a::zpz, 3 },{ "inc", &a::inc, &a::zpz, 5 },{ "???", &a::err, &a::imp, 5 },{ "inx", &a::inx, &a::imp, 2 },{ "sbc", &a::sbc, &a::imm, 2 },{ "nop", &a::nop, &a::imp, 2 },{ "???", &a::sbc, &a::imp, 2 },{ "cpx", &a::cpx, &a::abs, 4 },{ "sbc", &a::sbc, &a::abs, 4 },{ "inc", &a::inc, &a::abs, 6 },{ "???", &a::err, &a::imp, 6 },
        { "beq", &a::beq, &a::rel, 2 },{ "sbc", &a::sbc, &a::izy, 5 },{ "???", &a::err, &a::imp, 2 },{ "???", &a::err, &a::imp, 8 },{ "???", &a::nop, &a::imp, 4 },{ "sbc", &a::sbc, &a::zpx, 4 },{ "inc", &a::inc, &a::zpx, 6 },{ "???", &a::err, &a::imp, 6 },{ "sed", &a::sed, &a::imp, 2 },{ "sbc", &a::sbc, &a::aby, 4 },{ "nop", &a::nop, &a::imp, 2 },{ "???", &a::err, &a::imp, 7 },{ "???", &a::nop, &a::imp, 4 },{ "sbc", &a::sbc, &a::abx, 4 },{ "inc", &a::inc, &a::abx, 7 },{ "???", &a::err, &a::imp, 7 },
    };
} // used the method shown in olc's 6502 emulator video
cpu::~cpu(){}

u8 cpu::read(u16 address){
    return b->read(address);
}

void cpu::write(u16 address, u8 value){
    b->write(address, value);
}

u8 cpu::fetch(){
    if(lookup[opcode].addrMode != &cpu::imp)
        data = read(next);
    return data;
}

void cpu::reset(){
    u16 lo = read(0xFFFC);
    u16 hi = read(0xFFFD);

    pc = (hi << 8) | lo;

    a = x = y = 0;
    data = 0x00;
    next = nxtr = 0x0000;

    sp = 0xFD;
    st = 0b00100000;

    cycles = 8;
}

void cpu::irq(){
    if(!get(I)){

        write(0x0100 + sp--, (pc >> 8) & 0x00FF);
        write(0x0100 + sp--, pc & 0x00FF);

        assign(B, 0);
        assign(X, 1);
        assign(I, 1);

        write(0x0100 + sp--, st);

        pc = read(0xFFFE) << 8 | read(0xFFFF);

        cycles = 7;
    }
}

void cpu::nmi(){
    write(0x0100 + sp--, (pc >> 8) & 0x00FF);
    write(0x0100 + sp--, pc & 0x00FF);

    assign(B, 0);
    assign(X, 1);
    assign(I, 1);

    write(0x0100 + sp--, st);

    pc = read(0xFFFA) << 8 | read(0xFFFB);

    cycles = 7;
}

void cpu::clock(){
    if(cycles == 0){
        opcode = read(pc);
        pc++;
        assign(X,1);
    }
}


u8 cpu::imp(){
    data = a;
    return 0;
}

u8 cpu::imm(){
    next = pc++;
    return 0;
}

u8 cpu::zpz(){
    next = read(pc++) & 0x00FF;
    return 0;
}

u8 cpu::zpx(){
    next = (read(pc++)+x) & 0x00FF;
    return 0;
}

u8 cpu::zpy(){
    next = (read(pc++)+x) & 0x00FF;
    return 0;
}

u8 cpu::rel(){
    nxtr = read(pc++);
    if(nxtr & 0b10000000) nxtr |= 0xFF00; // last digit is 1 ie is negative
    return 0;
}

u8 cpu::abs(){
    next = read(pc++) | (read(pc++) << 8);
    return 0;
}

u8 cpu::abx(){
    u16 lo = read(pc++);
    u16 hi = read(pc++) << 8;
    next = hi | lo;
    next += x;

    if((next & 0xFF00) != hi) return 1; // checks page change
    return 0;
}

u8 cpu::aby(){
    u16 lo = read(pc++);
    u16 hi = read(pc++) << 8;
    next = hi | lo;
    next += y;

    if((next & 0xFF00) != hi) return 1; // checks page change
    return 0;
}

u8 cpu::ind(){
    u16 plo = read(pc++);
    u16 phi = read(pc++) << 8;
    u16 ptr = phi | plo;

    if(plo == 0x00FF) next = (read(phi) << 8) | read(ptr); // here, 0x00FF + 1 = 0x0000 thus, ptr + 1 = phi + plo + 1 = phi + 0x00FF + 1 = phi
    else    next = (read(ptr + 1) << 8) | read(ptr);

    return 0;
}

u8 cpu::izx(){
    u16 t = read(pc++);

    u16 lo = read((u16)((t + (u16)x) & 0x00FF));
    u16 hi = read((u16)((t + 1 + (u16)x) & 0x00FF)) << 8;

    next = hi | lo;
    return 0;
}

u8 cpu::izy(){
    u16 t = read(pc++);

    u16 lo = read(t & 0x00FF);
    u16 hi = read((t + 1) & 0x00FF) << 8;

    next = hi | lo;
    next += y;

    if((next & 0xFF00) != hi) return 1; // checks page change
    return 0;
}


u8 cpu::adc(){
	fetch();
	t = (u16)a + (u16)data + (u16)get(C);

	assign(C, t > 0x00FF);
	assign(Z, (t & 0x00FF) == 0);
	assign(V, (~((u16)a ^ (u16)data) & ((u16)a ^ (u16)t)) & 0x0080);
	assign(N, t & 0x80);
	a = t & 0x00FF;
	return 1;
}


u8 cpu::sbc(){
	fetch();
	
	u16 val = ((u16)data) ^ 0x00FF;
	t = (u16)a + val + (u16)get(C);

	assign(C, t & 0xFF00);
	assign(Z, ((t & 0x00FF) == 0));
	assign(V, (t ^ (u16)a) & (t ^ val) & 0x0080);
	assign(N, t & 0x0080);
	a = t & 0x00FF;
	return 1;
}

u8 cpu::nnd(){
	fetch();
	a = a & data;
	assign(Z, a == 0x00);
	assign(N, a & 0x80);
	return 1;
}

u8 cpu::asl(){
	fetch();
	t = (u16)data << 1;
	assign(C, (t & 0xFF00) > 0);
	assign(Z, (t & 0x00FF) == 0x00);
	assign(N, t & 0x80);
	if (lookup[opcode].addrMode == &cpu::imp)
		a = t & 0x00FF;
	else
		write(next, t & 0x00FF);
	return 0;
}

u8 cpu::bcc(){
	if (get(C) == 0)
	{
		cycles++;
		next = pc + nxtr;
		
		if((next & 0xFF00) != (pc & 0xFF00))
			cycles++;
		
		pc = next;
	}
	return 0;
}


u8 cpu::bcs(){
	if (get(C) == 1)
	{
		cycles++;
		next = pc + nxtr;

		if ((next & 0xFF00) != (pc & 0xFF00))
			cycles++;

		pc = next;
	}
	return 0;
}

u8 cpu::beq(){
	if (get(Z) == 1)
	{
		cycles++;
		next = pc + nxtr;

		if ((next & 0xFF00) != (pc & 0xFF00))
			cycles++;

		pc = next;
	}
	return 0;
}

u8 cpu::bit(){
	fetch();
	t = a & data;
	assign(Z, (t & 0x00FF) == 0x00);
	assign(N, data & (1 << 7));
	assign(V, data & (1 << 6));
	return 0;
}


u8 cpu::bmi(){
	if (get(N) == 1)
	{
		cycles++;
		next = pc + nxtr;

		if ((next & 0xFF00) != (pc & 0xFF00))
			cycles++;

		pc = next;
	}
	return 0;
}

u8 cpu::bne(){
	if (get(Z) == 0)
	{
		cycles++;
		next = pc + nxtr;

		if ((next & 0xFF00) != (pc & 0xFF00))
			cycles++;

		pc = next;
	}
	return 0;
}


u8 cpu::bpl(){
	if (get(N) == 0)
	{
		cycles++;
		next = pc + nxtr;

		if ((next & 0xFF00) != (pc & 0xFF00))
			cycles++;

		pc = next;
	}
	return 0;
}


u8 cpu::brk(){
	pc++;
	
	assign(I, 1);
	write(0x0100 + sp--, (pc >> 8) & 0x00FF);
	write(0x0100 + sp--, pc & 0x00FF);

	assign(B, 1);
	write(0x0100 + sp--, st);
	assign(B, 0);

	pc = (u16)read(0xFFFE) | ((u16)read(0xFFFF) << 8);
	return 0;
}

u8 cpu::bvc(){
	if (get(V) == 0){
		cycles++;
		next = pc + nxtr;

		if ((next & 0xFF00) != (pc & 0xFF00))
			cycles++;

		pc = next;
	}
	return 0;
}


u8 cpu::bvs(){
	if (get(V) == 1){
		cycles++;
		next = pc + nxtr;

		if ((next & 0xFF00) != (pc & 0xFF00))
			cycles++;

		pc = next;
	}
	return 0;
}


u8 cpu::clc(){
	assign(C, false);
	return 0;
}

u8 cpu::cld(){
	assign(D, false);
	return 0;
}

u8 cpu::cli(){
	assign(I, false);
	return 0;
}

u8 cpu::clv(){
	assign(V, false);
	return 0;
}

u8 cpu::cmp(){
	fetch();
	t = (u16)a - (u16)data;
	assign(C, a >= data);
	assign(Z, (t & 0x00FF) == 0x0000);
	assign(N, t & 0x0080);
	return 1;
}

u8 cpu::cpx(){
	fetch();
	t = (u16)x - (u16)data;
	assign(C, x >= data);
	assign(Z, (t & 0x00FF) == 0x0000);
	assign(N, t & 0x0080);
	return 0;
}

u8 cpu::cpy(){
	fetch();
	t = (u16)y - (u16)data;
	assign(C, y >= data);
	assign(Z, (t & 0x00FF) == 0x0000);
	assign(N, t & 0x0080);
	return 0;
}

u8 cpu::dec(){
	fetch();
	t = data - 1;
	write(next, t & 0x00FF);
	assign(Z, (t & 0x00FF) == 0x0000);
	assign(N, t & 0x0080);
	return 0;
}

u8 cpu::dex(){
	x--;
	assign(Z, x == 0x00);
	assign(N, x & 0x80);
	return 0;
}

u8 cpu::dey(){
	y--;
	assign(Z, y == 0x00);
	assign(N, y & 0x80);
	return 0;
}

u8 cpu::eor(){
	fetch();
	a = a ^ data;	
	assign(Z, a == 0x00);
	assign(N, a & 0x80);
	return 1;
}

u8 cpu::inc(){
	fetch();
	t = data + 1;
	write(next, t & 0x00FF);
	assign(Z, (t & 0x00FF) == 0x0000);
	assign(N, t & 0x0080);
	return 0;
}

u8 cpu::inx(){
	x++;
	assign(Z, x == 0x00);
	assign(N, x & 0x80);
	return 0;
}

u8 cpu::iny(){
	y++;
	assign(Z, y == 0x00);
	assign(N, y & 0x80);
	return 0;
}

u8 cpu::jmp(){
	pc = next;
	return 0;
}

u8 cpu::jsr(){
	pc--;

	write(0x0100 + sp--, (pc >> 8) & 0x00FF);
	write(0x0100 + sp--, pc & 0x00FF);

	pc = next;
	return 0;
}

u8 cpu::lda(){
	fetch();
	a = data;
	assign(Z, a == 0x00);
	assign(N, a & 0x80);
	return 1;
}

u8 cpu::ldx(){
	fetch();
	x = data;
	assign(Z, x == 0x00);
	assign(N, x & 0x80);
	return 1;
}

u8 cpu::ldy(){
	fetch();
	y = data;
	assign(Z, y == 0x00);
	assign(N, y & 0x80);
	return 1;
}

u8 cpu::lsr(){
	fetch();
	assign(C, data & 0x0001);
	t = data >> 1;	
	assign(Z, (t & 0x00FF) == 0x0000);
	assign(N, t & 0x0080);
	
	if (lookup[opcode].addrMode == &cpu::imp) a = t & 0x00FF;
	else write(next, t & 0x00FF);
	
	return 0;
}

u8 cpu::nop(){
	switch (opcode){
	case 0x1C:
	case 0x3C:
	case 0x5C:
	case 0x7C:
	case 0xDC:
	case 0xFC: return 1; break;
	}
	return 0;
}

u8 cpu::ora(){
	fetch();
	a = a | data;
	assign(Z, a == 0x00);
	assign(N, a & 0x80);
	return 1;
}

u8 cpu::pha(){
	write(0x0100 + sp--, a);
	return 0;
}

u8 cpu::php(){
	write(0x0100 + sp--, st | B | X);
	assign(B, 0);
	assign(X, 0);
	return 0;
}

u8 cpu::pla(){
	sp++;
	a = read(0x0100 + sp);
	assign(Z, a == 0x00);
	assign(N, a & 0x80);
	return 0;
}

u8 cpu::plp(){
	sp++;
	st = read(0x0100 + sp);
	assign(X, 1);
	return 0;
}

u8 cpu::rol(){
	fetch();
	t = (u16)(data << 1) | get(C);
	assign(C, t & 0xFF00);
	assign(Z, (t & 0x00FF) == 0x0000);
	assign(N, t & 0x0080);
	if (lookup[opcode].addrMode == &cpu::imp) a = t & 0x00FF;
	else write(next, t & 0x00FF);
	return 0;
}

u8 cpu::ror(){
	fetch();
	t = (u16)(get(C) << 7) | (data >> 1);
	assign(C, data & 0x01);
	assign(Z, (t & 0x00FF) == 0x00);
	assign(N, t & 0x0080);
	if (lookup[opcode].addrMode == &cpu::imp) a = t & 0x00FF;
	else write(next, t & 0x00FF);
	return 0;
}

u8 cpu::rti(){
	sp++;
	st = read(0x0100 + sp++);
	st &= ~B;
	st &= ~X;

	pc = (u16)read(0x0100 + sp++);
	pc |= (u16)read(0x0100 + sp) << 8;
	return 0;
}

u8 cpu::rts(){
	pc = (u16)read(0x0100 + ++sp);
	pc |= (u16)read(0x0100 + ++sp) << 8;
	pc++;
	return 0;
}

u8 cpu::sec(){
	assign(C, true);
	return 0;
}

u8 cpu::sed(){
	assign(D, true);
	return 0;
}
u8 cpu::sei(){
	assign(I, true);
	return 0;
}
u8 cpu::sta(){
	write(next, a);
	return 0;
}

u8 cpu::stx(){
	write(next, x);
	return 0;
}
u8 cpu::sty(){
	write(next, y);
	return 0;
}

u8 cpu::tax(){
	x = a;
	assign(Z, x == 0x00);
	assign(N, x & 0x80);
	return 0;
}

u8 cpu::tay(){
	y = a;
	assign(Z, y == 0x00);
	assign(N, y & 0x80);
	return 0;
}
u8 cpu::tsx(){
	x = sp;
	assign(Z, x == 0x00);
	assign(N, x & 0x80);
	return 0;
}
u8 cpu::txa(){
	a = x;
	assign(Z, a == 0x00);
	assign(N, a & 0x80);
	return 0;
}

u8 cpu::txs(){
	sp = x;
	return 0;
}

u8 cpu::tya(){
	a = y;
	assign(Z, a == 0x00);
	assign(N, a & 0x80);
	return 0;
}

u8 cpu::err(){
	return 0;
}
