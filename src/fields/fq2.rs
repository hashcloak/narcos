use ark_ff::{
    fields::fp2::{F2p, Fp2Config},
    Field, MontFp,
};

use crate::Fq;

pub type Fq2 = Fp2<Fp2Config>;

pub struct Fq2Config;

impl Fp2Config for Fq2Config {
    type Fp = Fq;

    /// NONRESIDUE = 7
    /// I.e. x^2 - 7 is irreducible in Fq
    const NONRESIDUE: Fq = MontFp!("7");

    const FROBENIUS_COEFF_FP2_C1: &'static [Fq] = &[
        // Fq(-1)**(((q^0) - 1) / 2)
        Fq::ONE,
        // Fq(-1)**(((q^1) - 1) / 2)
        MontFp!("458895104874816770381726682438519142702267188344382758289638358715359672478465538019994971562983093075690567609118224292550817460581580405075089814532843556636"),
    ];
}