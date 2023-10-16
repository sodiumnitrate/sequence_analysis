#include "include/pair_score.hpp"

PairScore::PairScore(std::string name_)
{
    name = name_;
    std::cout << name << std::endl;
    if (name.compare("blosum50") == 0) init_blosum50();
    else if (name.compare("blosum62") == 0) init_blosum62();
    else if (name.compare("blastn") == 0) init_blastn();
    else if (name.compare("levenshtein") == 0) return;
    else throw "unrecognized name";
}

void PairScore::init_blastn(){
    scores.clear();
    std::string letters = "ACGTUN";

    for(auto let1 : letters){
        for (auto let2 : letters){
            std::string pair;
            pair.push_back(let1);
            pair.push_back(let2);
            if (let1 == let2){
                scores[pair] = 1;
            }
            else{
                if (let1 == 'N' || let2 == 'N') scores[pair] = 1;
                else
                {
                    scores[pair] = -2;
                }
            }
        }
    }
}

void PairScore::init_blosum50(){
    scores.clear();
    // load data and populate map
    scores = {{"AA", 5},{"AR", -2},{"AN", -1},{"AD", -2},{"AC", -1},{"AQ", -1},{"AE", -1},{"AG", 0},{"AH", -2},{"AI", -1},{"AL", -2},{"AK", -1},{"AM", -1},{"AF", -3},{"AP", -1},{"AS", 1},{"AT", 0},{"AW", -3},{"AY", -2},{"AV", 0},{"AB", -2},{"AJ", -2},{"AZ", -1},{"AX", -1},{"A*", -5},
{"RA", -2},{"RR", 7},{"RN", -1},{"RD", -2},{"RC", -4},{"RQ", 1},{"RE", 0},{"RG", -3},{"RH", 0},{"RI", -4},{"RL", -3},{"RK", 3},{"RM", -2},{"RF", -3},{"RP", -3},{"RS", -1},{"RT", -1},{"RW", -3},{"RY", -1},{"RV", -3},{"RB", -1},{"RJ", -3},{"RZ", 0},{"RX", -1},{"R*", -5},
{"NA", -1},{"NR", -1},{"NN", 7},{"ND", 2},{"NC", -2},{"NQ", 0},{"NE", 0},{"NG", 0},{"NH", 1},{"NI", -3},{"NL", -4},{"NK", 0},{"NM", -2},{"NF", -4},{"NP", -2},{"NS", 1},{"NT", 0},{"NW", -4},{"NY", -2},{"NV", -3},{"NB", 5},{"NJ", -4},{"NZ", 0},{"NX", -1},{"N*", -5},
{"DA", -2},{"DR", -2},{"DN", 2},{"DD", 8},{"DC", -4},{"DQ", 0},{"DE", 2},{"DG", -1},{"DH", -1},{"DI", -4},{"DL", -4},{"DK", -1},{"DM", -4},{"DF", -5},{"DP", -1},{"DS", 0},{"DT", -1},{"DW", -5},{"DY", -3},{"DV", -4},{"DB", 6},{"DJ", -4},{"DZ", 1},{"DX", -1},{"D*", -5},
{"CA", -1},{"CR", -4},{"CN", -2},{"CD", -4},{"CC", 13},{"CQ", -3},{"CE", -3},{"CG", -3},{"CH", -3},{"CI", -2},{"CL", -2},{"CK", -3},{"CM", -2},{"CF", -2},{"CP", -4},{"CS", -1},{"CT", -1},{"CW", -5},{"CY", -3},{"CV", -1},{"CB", -3},{"CJ", -2},{"CZ", -3},{"CX", -1},{"C*", -5},
{"QA", -1},{"QR", 1},{"QN", 0},{"QD", 0},{"QC", -3},{"QQ", 7},{"QE", 2},{"QG", -2},{"QH", 1},{"QI", -3},{"QL", -2},{"QK", 2},{"QM", 0},{"QF", -4},{"QP", -1},{"QS", 0},{"QT", -1},{"QW", -1},{"QY", -1},{"QV", -3},{"QB", 0},{"QJ", -3},{"QZ", 4},{"QX", -1},{"Q*", -5},
{"EA", -1},{"ER", 0},{"EN", 0},{"ED", 2},{"EC", -3},{"EQ", 2},{"EE", 6},{"EG", -3},{"EH", 0},{"EI", -4},{"EL", -3},{"EK", 1},{"EM", -2},{"EF", -3},{"EP", -1},{"ES", -1},{"ET", -1},{"EW", -3},{"EY", -2},{"EV", -3},{"EB", 1},{"EJ", -3},{"EZ", 5},{"EX", -1},{"E*", -5},
{"GA", 0},{"GR", -3},{"GN", 0},{"GD", -1},{"GC", -3},{"GQ", -2},{"GE", -3},{"GG", 8},{"GH", -2},{"GI", -4},{"GL", -4},{"GK", -2},{"GM", -3},{"GF", -4},{"GP", -2},{"GS", 0},{"GT", -2},{"GW", -3},{"GY", -3},{"GV", -4},{"GB", -1},{"GJ", -4},{"GZ", -2},{"GX", -1},{"G*", -5},
{"HA", -2},{"HR", 0},{"HN", 1},{"HD", -1},{"HC", -3},{"HQ", 1},{"HE", 0},{"HG", -2},{"HH", 10},{"HI", -4},{"HL", -3},{"HK", 0},{"HM", -1},{"HF", -1},{"HP", -2},{"HS", -1},{"HT", -2},{"HW", -3},{"HY", 2},{"HV", -4},{"HB", 0},{"HJ", -3},{"HZ", 0},{"HX", -1},{"H*", -5},
{"IA", -1},{"IR", -4},{"IN", -3},{"ID", -4},{"IC", -2},{"IQ", -3},{"IE", -4},{"IG", -4},{"IH", -4},{"II", 5},{"IL", 2},{"IK", -3},{"IM", 2},{"IF", 0},{"IP", -3},{"IS", -3},{"IT", -1},{"IW", -3},{"IY", -1},{"IV", 4},{"IB", -4},{"IJ", 4},{"IZ", -3},{"IX", -1},{"I*", -5},
{"LA", -2},{"LR", -3},{"LN", -4},{"LD", -4},{"LC", -2},{"LQ", -2},{"LE", -3},{"LG", -4},{"LH", -3},{"LI", 2},{"LL", 5},{"LK", -3},{"LM", 3},{"LF", 1},{"LP", -4},{"LS", -3},{"LT", -1},{"LW", -2},{"LY", -1},{"LV", 1},{"LB", -4},{"LJ", 4},{"LZ", -3},{"LX", -1},{"L*", -5},
{"KA", -1},{"KR", 3},{"KN", 0},{"KD", -1},{"KC", -3},{"KQ", 2},{"KE", 1},{"KG", -2},{"KH", 0},{"KI", -3},{"KL", -3},{"KK", 6},{"KM", -2},{"KF", -4},{"KP", -1},{"KS", 0},{"KT", -1},{"KW", -3},{"KY", -2},{"KV", -3},{"KB", 0},{"KJ", -3},{"KZ", 1},{"KX", -1},{"K*", -5},
{"MA", -1},{"MR", -2},{"MN", -2},{"MD", -4},{"MC", -2},{"MQ", 0},{"ME", -2},{"MG", -3},{"MH", -1},{"MI", 2},{"ML", 3},{"MK", -2},{"MM", 7},{"MF", 0},{"MP", -3},{"MS", -2},{"MT", -1},{"MW", -1},{"MY", 0},{"MV", 1},{"MB", -3},{"MJ", 2},{"MZ", -1},{"MX", -1},{"M*", -5},
{"FA", -3},{"FR", -3},{"FN", -4},{"FD", -5},{"FC", -2},{"FQ", -4},{"FE", -3},{"FG", -4},{"FH", -1},{"FI", 0},{"FL", 1},{"FK", -4},{"FM", 0},{"FF", 8},{"FP", -4},{"FS", -3},{"FT", -2},{"FW", 1},{"FY", 4},{"FV", -1},{"FB", -4},{"FJ", 1},{"FZ", -4},{"FX", -1},{"F*", -5},
{"PA", -1},{"PR", -3},{"PN", -2},{"PD", -1},{"PC", -4},{"PQ", -1},{"PE", -1},{"PG", -2},{"PH", -2},{"PI", -3},{"PL", -4},{"PK", -1},{"PM", -3},{"PF", -4},{"PP", 10},{"PS", -1},{"PT", -1},{"PW", -4},{"PY", -3},{"PV", -3},{"PB", -2},{"PJ", -3},{"PZ", -1},{"PX", -1},{"P*", -5},
{"SA", 1},{"SR", -1},{"SN", 1},{"SD", 0},{"SC", -1},{"SQ", 0},{"SE", -1},{"SG", 0},{"SH", -1},{"SI", -3},{"SL", -3},{"SK", 0},{"SM", -2},{"SF", -3},{"SP", -1},{"SS", 5},{"ST", 2},{"SW", -4},{"SY", -2},{"SV", -2},{"SB", 0},{"SJ", -3},{"SZ", 0},{"SX", -1},{"S*", -5},
{"TA", 0},{"TR", -1},{"TN", 0},{"TD", -1},{"TC", -1},{"TQ", -1},{"TE", -1},{"TG", -2},{"TH", -2},{"TI", -1},{"TL", -1},{"TK", -1},{"TM", -1},{"TF", -2},{"TP", -1},{"TS", 2},{"TT", 5},{"TW", -3},{"TY", -2},{"TV", 0},{"TB", 0},{"TJ", -1},{"TZ", -1},{"TX", -1},{"T*", -5},
{"WA", -3},{"WR", -3},{"WN", -4},{"WD", -5},{"WC", -5},{"WQ", -1},{"WE", -3},{"WG", -3},{"WH", -3},{"WI", -3},{"WL", -2},{"WK", -3},{"WM", -1},{"WF", 1},{"WP", -4},{"WS", -4},{"WT", -3},{"WW", 15},{"WY", 2},{"WV", -3},{"WB", -5},{"WJ", -2},{"WZ", -2},{"WX", -1},{"W*", -5},
{"YA", -2},{"YR", -1},{"YN", -2},{"YD", -3},{"YC", -3},{"YQ", -1},{"YE", -2},{"YG", -3},{"YH", 2},{"YI", -1},{"YL", -1},{"YK", -2},{"YM", 0},{"YF", 4},{"YP", -3},{"YS", -2},{"YT", -2},{"YW", 2},{"YY", 8},{"YV", -1},{"YB", -3},{"YJ", -1},{"YZ", -2},{"YX", -1},{"Y*", -5},
{"VA", 0},{"VR", -3},{"VN", -3},{"VD", -4},{"VC", -1},{"VQ", -3},{"VE", -3},{"VG", -4},{"VH", -4},{"VI", 4},{"VL", 1},{"VK", -3},{"VM", 1},{"VF", -1},{"VP", -3},{"VS", -2},{"VT", 0},{"VW", -3},{"VY", -1},{"VV", 5},{"VB", -3},{"VJ", 2},{"VZ", -3},{"VX", -1},{"V*", -5},
{"BA", -2},{"BR", -1},{"BN", 5},{"BD", 6},{"BC", -3},{"BQ", 0},{"BE", 1},{"BG", -1},{"BH", 0},{"BI", -4},{"BL", -4},{"BK", 0},{"BM", -3},{"BF", -4},{"BP", -2},{"BS", 0},{"BT", 0},{"BW", -5},{"BY", -3},{"BV", -3},{"BB", 6},{"BJ", -4},{"BZ", 1},{"BX", -1},{"B*", -5},
{"JA", -2},{"JR", -3},{"JN", -4},{"JD", -4},{"JC", -2},{"JQ", -3},{"JE", -3},{"JG", -4},{"JH", -3},{"JI", 4},{"JL", 4},{"JK", -3},{"JM", 2},{"JF", 1},{"JP", -3},{"JS", -3},{"JT", -1},{"JW", -2},{"JY", -1},{"JV", 2},{"JB", -4},{"JJ", 4},{"JZ", -3},{"JX", -1},{"J*", -5},
{"ZA", -1},{"ZR", 0},{"ZN", 0},{"ZD", 1},{"ZC", -3},{"ZQ", 4},{"ZE", 5},{"ZG", -2},{"ZH", 0},{"ZI", -3},{"ZL", -3},{"ZK", 1},{"ZM", -1},{"ZF", -4},{"ZP", -1},{"ZS", 0},{"ZT", -1},{"ZW", -2},{"ZY", -2},{"ZV", -3},{"ZB", 1},{"ZJ", -3},{"ZZ", 5},{"ZX", -1},{"Z*", -5},
{"XA", -1},{"XR", -1},{"XN", -1},{"XD", -1},{"XC", -1},{"XQ", -1},{"XE", -1},{"XG", -1},{"XH", -1},{"XI", -1},{"XL", -1},{"XK", -1},{"XM", -1},{"XF", -1},{"XP", -1},{"XS", -1},{"XT", -1},{"XW", -1},{"XY", -1},{"XV", -1},{"XB", -1},{"XJ", -1},{"XZ", -1},{"XX", -1},{"X*", -5},
{"*A", -5},{"*R", -5},{"*N", -5},{"*D", -5},{"*C", -5},{"*Q", -5},{"*E", -5},{"*G", -5},{"*H", -5},{"*I", -5},{"*L", -5},{"*K", -5},{"*M", -5},{"*F", -5},{"*P", -5},{"*S", -5},{"*T", -5},{"*W", -5},{"*Y", -5},{"*V", -5},{"*B", -5},{"*J", -5},{"*Z", -5},{"*X", -5},{"**", 1}};
}

void PairScore::init_blosum62(){
    scores.clear();
    scores = {
        {"AA", 4},{"AR", -1},{"AN", -2},{"AD", -2},{"AC", 0},{"AQ", -1},{"AE", -1},{"AG", 0},{"AH", -2},{"AI", -1},{"AL", -1},{"AK", -1},{"AM", -1},{"AF", -2},{"AP", -1},{"AS", 1},{"AT", 0},{"AW", -3},{"AY", -2},{"AV", 0},{"AB", -2},{"AJ", -1},{"AZ", -1},{"AX", -1},{"A*", -4},
{"RA", -1},{"RR", 5},{"RN", 0},{"RD", -2},{"RC", -3},{"RQ", 1},{"RE", 0},{"RG", -2},{"RH", 0},{"RI", -3},{"RL", -2},{"RK", 2},{"RM", -1},{"RF", -3},{"RP", -2},{"RS", -1},{"RT", -1},{"RW", -3},{"RY", -2},{"RV", -3},{"RB", -1},{"RJ", -2},{"RZ", 0},{"RX", -1},{"R*", -4},
{"NA", -2},{"NR", 0},{"NN", 6},{"ND", 1},{"NC", -3},{"NQ", 0},{"NE", 0},{"NG", 0},{"NH", 1},{"NI", -3},{"NL", -3},{"NK", 0},{"NM", -2},{"NF", -3},{"NP", -2},{"NS", 1},{"NT", 0},{"NW", -4},{"NY", -2},{"NV", -3},{"NB", 4},{"NJ", -3},{"NZ", 0},{"NX", -1},{"N*", -4},
{"DA", -2},{"DR", -2},{"DN", 1},{"DD", 6},{"DC", -3},{"DQ", 0},{"DE", 2},{"DG", -1},{"DH", -1},{"DI", -3},{"DL", -4},{"DK", -1},{"DM", -3},{"DF", -3},{"DP", -1},{"DS", 0},{"DT", -1},{"DW", -4},{"DY", -3},{"DV", -3},{"DB", 4},{"DJ", -3},{"DZ", 1},{"DX", -1},{"D*", -4},
{"CA", 0},{"CR", -3},{"CN", -3},{"CD", -3},{"CC", 9},{"CQ", -3},{"CE", -4},{"CG", -3},{"CH", -3},{"CI", -1},{"CL", -1},{"CK", -3},{"CM", -1},{"CF", -2},{"CP", -3},{"CS", -1},{"CT", -1},{"CW", -2},{"CY", -2},{"CV", -1},{"CB", -3},{"CJ", -1},{"CZ", -3},{"CX", -1},{"C*", -4},
{"QA", -1},{"QR", 1},{"QN", 0},{"QD", 0},{"QC", -3},{"QQ", 5},{"QE", 2},{"QG", -2},{"QH", 0},{"QI", -3},{"QL", -2},{"QK", 1},{"QM", 0},{"QF", -3},{"QP", -1},{"QS", 0},{"QT", -1},{"QW", -2},{"QY", -1},{"QV", -2},{"QB", 0},{"QJ", -2},{"QZ", 4},{"QX", -1},{"Q*", -4},
{"EA", -1},{"ER", 0},{"EN", 0},{"ED", 2},{"EC", -4},{"EQ", 2},{"EE", 5},{"EG", -2},{"EH", 0},{"EI", -3},{"EL", -3},{"EK", 1},{"EM", -2},{"EF", -3},{"EP", -1},{"ES", 0},{"ET", -1},{"EW", -3},{"EY", -2},{"EV", -2},{"EB", 1},{"EJ", -3},{"EZ", 4},{"EX", -1},{"E*", -4},
{"GA", 0},{"GR", -2},{"GN", 0},{"GD", -1},{"GC", -3},{"GQ", -2},{"GE", -2},{"GG", 6},{"GH", -2},{"GI", -4},{"GL", -4},{"GK", -2},{"GM", -3},{"GF", -3},{"GP", -2},{"GS", 0},{"GT", -2},{"GW", -2},{"GY", -3},{"GV", -3},{"GB", -1},{"GJ", -4},{"GZ", -2},{"GX", -1},{"G*", -4},
{"HA", -2},{"HR", 0},{"HN", 1},{"HD", -1},{"HC", -3},{"HQ", 0},{"HE", 0},{"HG", -2},{"HH", 8},{"HI", -3},{"HL", -3},{"HK", -1},{"HM", -2},{"HF", -1},{"HP", -2},{"HS", -1},{"HT", -2},{"HW", -2},{"HY", 2},{"HV", -3},{"HB", 0},{"HJ", -3},{"HZ", 0},{"HX", -1},{"H*", -4},
{"IA", -1},{"IR", -3},{"IN", -3},{"ID", -3},{"IC", -1},{"IQ", -3},{"IE", -3},{"IG", -4},{"IH", -3},{"II", 4},{"IL", 2},{"IK", -3},{"IM", 1},{"IF", 0},{"IP", -3},{"IS", -2},{"IT", -1},{"IW", -3},{"IY", -1},{"IV", 3},{"IB", -3},{"IJ", 3},{"IZ", -3},{"IX", -1},{"I*", -4},
{"LA", -1},{"LR", -2},{"LN", -3},{"LD", -4},{"LC", -1},{"LQ", -2},{"LE", -3},{"LG", -4},{"LH", -3},{"LI", 2},{"LL", 4},{"LK", -2},{"LM", 2},{"LF", 0},{"LP", -3},{"LS", -2},{"LT", -1},{"LW", -2},{"LY", -1},{"LV", 1},{"LB", -4},{"LJ", 3},{"LZ", -3},{"LX", -1},{"L*", -4},
{"KA", -1},{"KR", 2},{"KN", 0},{"KD", -1},{"KC", -3},{"KQ", 1},{"KE", 1},{"KG", -2},{"KH", -1},{"KI", -3},{"KL", -2},{"KK", 5},{"KM", -1},{"KF", -3},{"KP", -1},{"KS", 0},{"KT", -1},{"KW", -3},{"KY", -2},{"KV", -2},{"KB", 0},{"KJ", -3},{"KZ", 1},{"KX", -1},{"K*", -4},
{"MA", -1},{"MR", -1},{"MN", -2},{"MD", -3},{"MC", -1},{"MQ", 0},{"ME", -2},{"MG", -3},{"MH", -2},{"MI", 1},{"ML", 2},{"MK", -1},{"MM", 5},{"MF", 0},{"MP", -2},{"MS", -1},{"MT", -1},{"MW", -1},{"MY", -1},{"MV", 1},{"MB", -3},{"MJ", 2},{"MZ", -1},{"MX", -1},{"M*", -4},
{"FA", -2},{"FR", -3},{"FN", -3},{"FD", -3},{"FC", -2},{"FQ", -3},{"FE", -3},{"FG", -3},{"FH", -1},{"FI", 0},{"FL", 0},{"FK", -3},{"FM", 0},{"FF", 6},{"FP", -4},{"FS", -2},{"FT", -2},{"FW", 1},{"FY", 3},{"FV", -1},{"FB", -3},{"FJ", 0},{"FZ", -3},{"FX", -1},{"F*", -4},
{"PA", -1},{"PR", -2},{"PN", -2},{"PD", -1},{"PC", -3},{"PQ", -1},{"PE", -1},{"PG", -2},{"PH", -2},{"PI", -3},{"PL", -3},{"PK", -1},{"PM", -2},{"PF", -4},{"PP", 7},{"PS", -1},{"PT", -1},{"PW", -4},{"PY", -3},{"PV", -2},{"PB", -2},{"PJ", -3},{"PZ", -1},{"PX", -1},{"P*", -4},
{"SA", 1},{"SR", -1},{"SN", 1},{"SD", 0},{"SC", -1},{"SQ", 0},{"SE", 0},{"SG", 0},{"SH", -1},{"SI", -2},{"SL", -2},{"SK", 0},{"SM", -1},{"SF", -2},{"SP", -1},{"SS", 4},{"ST", 1},{"SW", -3},{"SY", -2},{"SV", -2},{"SB", 0},{"SJ", -2},{"SZ", 0},{"SX", -1},{"S*", -4},
{"TA", 0},{"TR", -1},{"TN", 0},{"TD", -1},{"TC", -1},{"TQ", -1},{"TE", -1},{"TG", -2},{"TH", -2},{"TI", -1},{"TL", -1},{"TK", -1},{"TM", -1},{"TF", -2},{"TP", -1},{"TS", 1},{"TT", 5},{"TW", -2},{"TY", -2},{"TV", 0},{"TB", -1},{"TJ", -1},{"TZ", -1},{"TX", -1},{"T*", -4},
{"WA", -3},{"WR", -3},{"WN", -4},{"WD", -4},{"WC", -2},{"WQ", -2},{"WE", -3},{"WG", -2},{"WH", -2},{"WI", -3},{"WL", -2},{"WK", -3},{"WM", -1},{"WF", 1},{"WP", -4},{"WS", -3},{"WT", -2},{"WW", 11},{"WY", 2},{"WV", -3},{"WB", -4},{"WJ", -2},{"WZ", -2},{"WX", -1},{"W*", -4},
{"YA", -2},{"YR", -2},{"YN", -2},{"YD", -3},{"YC", -2},{"YQ", -1},{"YE", -2},{"YG", -3},{"YH", 2},{"YI", -1},{"YL", -1},{"YK", -2},{"YM", -1},{"YF", 3},{"YP", -3},{"YS", -2},{"YT", -2},{"YW", 2},{"YY", 7},{"YV", -1},{"YB", -3},{"YJ", -1},{"YZ", -2},{"YX", -1},{"Y*", -4},
{"VA", 0},{"VR", -3},{"VN", -3},{"VD", -3},{"VC", -1},{"VQ", -2},{"VE", -2},{"VG", -3},{"VH", -3},{"VI", 3},{"VL", 1},{"VK", -2},{"VM", 1},{"VF", -1},{"VP", -2},{"VS", -2},{"VT", 0},{"VW", -3},{"VY", -1},{"VV", 4},{"VB", -3},{"VJ", 2},{"VZ", -2},{"VX", -1},{"V*", -4},
{"BA", -2},{"BR", -1},{"BN", 4},{"BD", 4},{"BC", -3},{"BQ", 0},{"BE", 1},{"BG", -1},{"BH", 0},{"BI", -3},{"BL", -4},{"BK", 0},{"BM", -3},{"BF", -3},{"BP", -2},{"BS", 0},{"BT", -1},{"BW", -4},{"BY", -3},{"BV", -3},{"BB", 4},{"BJ", -3},{"BZ", 0},{"BX", -1},{"B*", -4},
{"JA", -1},{"JR", -2},{"JN", -3},{"JD", -3},{"JC", -1},{"JQ", -2},{"JE", -3},{"JG", -4},{"JH", -3},{"JI", 3},{"JL", 3},{"JK", -3},{"JM", 2},{"JF", 0},{"JP", -3},{"JS", -2},{"JT", -1},{"JW", -2},{"JY", -1},{"JV", 2},{"JB", -3},{"JJ", 3},{"JZ", -3},{"JX", -1},{"J*", -4},
{"ZA", -1},{"ZR", 0},{"ZN", 0},{"ZD", 1},{"ZC", -3},{"ZQ", 4},{"ZE", 4},{"ZG", -2},{"ZH", 0},{"ZI", -3},{"ZL", -3},{"ZK", 1},{"ZM", -1},{"ZF", -3},{"ZP", -1},{"ZS", 0},{"ZT", -1},{"ZW", -2},{"ZY", -2},{"ZV", -2},{"ZB", 0},{"ZJ", -3},{"ZZ", 4},{"ZX", -1},{"Z*", -4},
{"XA", -1},{"XR", -1},{"XN", -1},{"XD", -1},{"XC", -1},{"XQ", -1},{"XE", -1},{"XG", -1},{"XH", -1},{"XI", -1},{"XL", -1},{"XK", -1},{"XM", -1},{"XF", -1},{"XP", -1},{"XS", -1},{"XT", -1},{"XW", -1},{"XY", -1},{"XV", -1},{"XB", -1},{"XJ", -1},{"XZ", -1},{"XX", -1},{"X*", -4},
{"*A", -4},{"*R", -4},{"*N", -4},{"*D", -4},{"*C", -4},{"*Q", -4},{"*E", -4},{"*G", -4},{"*H", -4},{"*I", -4},{"*L", -4},{"*K", -4},{"*M", -4},{"*F", -4},{"*P", -4},{"*S", -4},{"*T", -4},{"*W", -4},{"*Y", -4},{"*V", -4},{"*B", -4},{"*J", -4},{"*Z", -4},{"*X", -4},{"**", 1}
    };
}

int PairScore::query(char res1, char res2){
    std::string pair;
    pair.push_back(res1);
    pair.push_back(res2);
    if (scores.find(pair) == scores.end()) {
        std::cout << "char not found: " << pair << std::endl;
        throw "char not found";
    }
    return scores[pair];
}

