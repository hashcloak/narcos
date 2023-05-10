//! This module implements Narcos (CP6-529), a pairing-friendly curve that embeds the NIST standard P-256 curve
//! This curve has an estimated bits of security of 126 bits. For more details on how this was determined, see the candidates folder in this repository

mod curves;
mod fields;

pub use curves::*;
pub use fields::*;