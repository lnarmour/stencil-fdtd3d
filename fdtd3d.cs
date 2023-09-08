prog = ReadAlphabets("./fdtd3d.ab");
system = "fdtd3d";
 
outDir = "./codegen";

setSpaceTimeMap(prog, system, "ex", "(i,j,k->1,0,0,i,j,k)");
setSpaceTimeMap(prog, system, "ey", "(i,j,k->1,0,0,i,j,k)");
setSpaceTimeMap(prog, system, "ez", "(i,j,k->1,0,0,i,j,k)");

setSpaceTimeMap(prog, system, "Chxey", "(->2,0,0,0,0,0)");
setSpaceTimeMap(prog, system, "Chxez", "(->2,0,0,0,0,0)");
setSpaceTimeMap(prog, system, "Chyex", "(->2,0,0,0,0,0)");
setSpaceTimeMap(prog, system, "Chyez", "(->2,0,0,0,0,0)");
setSpaceTimeMap(prog, system, "Chzex", "(->2,0,0,0,0,0)");
setSpaceTimeMap(prog, system, "Chzey", "(->2,0,0,0,0,0)");
setSpaceTimeMap(prog, system, "Cexhy", "(i,j,k->3,0,0,i,j,k)");
setSpaceTimeMap(prog, system, "Cexhz", "(i,j,k->3,0,0,i,j,k)");
setSpaceTimeMap(prog, system, "Ceyhx", "(i,j,k->4,0,0,i,j,k)");
setSpaceTimeMap(prog, system, "Ceyhz", "(i,j,k->4,0,0,i,j,k)");
setSpaceTimeMap(prog, system, "Cezhx", "(i,j,k->5,0,0,i,j,k)");
setSpaceTimeMap(prog, system, "Cezhy", "(i,j,k->5,0,0,i,j,k)");

setSpaceTimeMap(prog, system, "Hx", "(t,i,j,k->6,t,0,i,j,k)");
setSpaceTimeMap(prog, system, "Hy", "(t,i,j,k->6,t,0,i,j,k)");
setSpaceTimeMap(prog, system, "Hz", "(t,i,j,k->6,t,0,i,j,k)");
setSpaceTimeMap(prog, system, "S",  "(t,i,j,k->6,t,1,i,j,k)");
setSpaceTimeMap(prog, system, "Ex", "(t,i,j,k->6,t,2,i,j,k)");
setSpaceTimeMap(prog, system, "Ey", "(t,i,j,k->6,t,2,i,j,k)");
setSpaceTimeMap(prog, system, "Ez", "(t,i,j,k->6,t,2,i,j,k)");
setSpaceTimeMap(prog, system, "V",  "(t,i,j,k->6,t,3,i,j,k)");


# MEMORY LAYOUTS

setMemorySpace(prog, system, "Cexhy", "Cexhy");
setMemoryMap(prog, system,   "Cexhy", "Cexhy", "(i,j,k->j)");
setMemorySpace(prog, system, "Cexhz", "Cexhz");
setMemoryMap(prog, system,   "Cexhz", "Cexhz", "(i,j,k->j)");
setMemorySpace(prog, system, "Ceyhx", "Ceyhx");
setMemoryMap(prog, system,   "Ceyhx", "Ceyhx", "(i,j,k->j)");
setMemorySpace(prog, system, "Ceyhz", "Ceyhz");
setMemoryMap(prog, system,   "Ceyhz", "Ceyhz", "(i,j,k->j)");
setMemorySpace(prog, system, "Cezhx", "Cezhx");
setMemoryMap(prog, system,   "Cezhx", "Cezhx", "(i,j,k->j)");
setMemorySpace(prog, system, "Cezhy", "Cezhy");
setMemoryMap(prog, system,   "Cezhy", "Cezhy", "(i,j,k->j)");

setMemorySpace(prog, system, "S", "S");
setMemoryMap(prog, system,   "S", "S", "(t,i,j,k->i,j,k)");

setMemorySpace(prog, system, "V", "V");
setMemoryMap(prog, system,   "V", "V", "(t,i,j,k->)");

setMemorySpace(prog, system, "Hx", "Hx");
setMemoryMap(prog, system,   "Hx", "Hx", "(t,i,j,k->t,i,j,k)", "2,");
setMemorySpace(prog, system, "Hy", "Hy");
setMemoryMap(prog, system,   "Hy", "Hy", "(t,i,j,k->t,i,j,k)", "2,");
setMemorySpace(prog, system, "Hz", "Hz");
setMemoryMap(prog, system,   "Hz", "Hz", "(t,i,j,k->t,i,j,k)", "2,");

setMemorySpace(prog, system, "Ex", "Ex");
setMemoryMap(prog, system,   "Ex", "Ex", "(t,i,j,k->t,i,j,k)", "2,");
setMemorySpace(prog, system, "Ey", "Ey");
setMemoryMap(prog, system,   "Ey", "Ey", "(t,i,j,k->t,i,j,k)", "2,");
setMemorySpace(prog, system, "Ez", "Ez");
setMemoryMap(prog, system,   "Ez", "Ez", "(t,i,j,k->t,i,j,k)", "2,");

#setParallel(prog, system, "", "3");

Show(prog);

generateScheduledCode(prog, system, outDir);
