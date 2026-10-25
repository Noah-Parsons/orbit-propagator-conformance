// Drives Basilisk's LambertSolver module directly (no SWIG): reads problems from stdin,
// one per line: method(0=Gooding,1=Izzo) N mu tof r1x r1y r1z r2x r2y r2z
// prints: valid validSol2 v1x v1y v1z v1Sol2x v1Sol2y v1Sol2z
#include "fswAlgorithms/orbitControl/lambertSolver/lambertSolver.h"
#include <cstdio>
#include <iostream>
int main() {
    int method, N; double mu, tof, a[6];
    while (std::cin >> method >> N >> mu >> tof >> a[0] >> a[1] >> a[2] >> a[3] >> a[4] >> a[5]) {
        LambertSolver mod;
        Message<LambertProblemMsgPayload> in;
        LambertProblemMsgPayload p = in.zeroMsgPayload;
        p.solverMethod = method == 0 ? GOODING : IZZO;
        for (int i = 0; i < 3; i++) { p.r1_N[i] = a[i]; p.r2_N[i] = a[3 + i]; }
        p.transferTime = tof; p.mu = mu; p.numRevolutions = N;
        in.write(&p, 0, 0);
        mod.lambertProblemInMsg.subscribeTo(&in);
        mod.Reset(0);
        mod.UpdateState(0);
        ReadFunctor<LambertSolutionMsgPayload> rd = mod.lambertSolutionOutMsg.addSubscriber();
        LambertSolutionMsgPayload s = rd();
        std::printf("%d %d %.15g %.15g %.15g %.15g %.15g %.15g\n", s.valid, s.validSol2, s.v1_N[0], s.v1_N[1], s.v1_N[2],
                    s.v1Sol2_N[0], s.v1Sol2_N[1], s.v1Sol2_N[2]);
    }
}
