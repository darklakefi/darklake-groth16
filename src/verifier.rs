use ark_ec::{pairing::Pairing, AffineRepr, CurveGroup};
use ark_ff::PrimeField;

use crate::{r1cs_to_qap::R1CSToQAP, Groth16};

use super::{PreparedVerifyingKey, Proof, VerifyingKey};

use ark_relations::r1cs::{Result as R1CSResult, SynthesisError};

use core::ops::{AddAssign, Neg};

/// Prepare the verifying key `vk` for use in proof verification.
pub fn prepare_verifying_key<E: Pairing>(vk: &VerifyingKey<E>) -> PreparedVerifyingKey<E> {
    PreparedVerifyingKey {
        vk: vk.clone(),
        alpha_g1_beta_g2: E::pairing(vk.alpha_g1, vk.beta_g2).0,
        gamma_g2_neg_pc: vk.gamma_g2.into_group().neg().into_affine().into(),
        delta_g2_neg_pc: vk.delta_g2.into_group().neg().into_affine().into(),
    }
}

impl<E: Pairing, QAP: R1CSToQAP> Groth16<E, QAP> {
    /// Prepare static and variable inputs for verification without modifying the VerifyingKey structure
    pub fn prepare_inputs_with_variables(
        pvk: &PreparedVerifyingKey<E>,
        static_inputs: &[E::ScalarField],
        variable_inputs: &[E::ScalarField],
    ) -> R1CSResult<E::G1> {
        // Calculate the expected total inputs length
        let expected_total_inputs = static_inputs.len() + variable_inputs.len();
        
        // Check that the number of inputs matches the gamma_abc_g1 size
        if expected_total_inputs + 1 != (pvk.vk.gamma_abc_g1_static.len() + pvk.vk.gamma_abc_g1_variable.len()) {
            return Err(SynthesisError::MalformedVerifyingKey);
        }

        // Start with the constant term
        let mut g_ic = pvk.vk.gamma_abc_g1_static[0].into_group();
        
        // Add static inputs first
        for (i, static_input) in static_inputs.iter().enumerate() {
            let idx = i + 1; // +1 for the constant term
            g_ic.add_assign(&pvk.vk.gamma_abc_g1_static[idx].mul_bigint(static_input.into_bigint()));
        }
        
        // Add variable inputs
        let static_offset = static_inputs.len() + 1; // +1 for the constant term
        for (i, variable_input) in variable_inputs.iter().enumerate() {
            let idx = static_offset + i;
            g_ic.add_assign(&pvk.vk.gamma_abc_g1_variable[idx].mul_bigint(variable_input.into_bigint()));
        }
        
        Ok(g_ic)
    }

    /// Verify a proof with separate static and variable inputs
    pub fn verify_with_variables(
        pvk: &PreparedVerifyingKey<E>,
        proof: &Proof<E>,
        static_inputs: &[E::ScalarField],
        variable_inputs: &[E::ScalarField],
    ) -> R1CSResult<bool> {
        // Prepare the input aggregation
        let prepared_inputs = Self::prepare_inputs_with_variables(pvk, static_inputs, variable_inputs)?;
        
        // Use standard verification with prepared inputs
        Self::verify_proof_with_prepared_inputs(pvk, proof, &prepared_inputs)
    }
        
    /// Verify a Groth16 proof `proof` against the prepared verification key `pvk` and prepared public
    /// inputs. This should be preferred over [`verify_proof`] if the instance's public inputs are
    /// known in advance.
    pub fn verify_proof_with_prepared_inputs(
        pvk: &PreparedVerifyingKey<E>,
        proof: &Proof<E>,
        prepared_inputs: &E::G1,
    ) -> R1CSResult<bool> {
        let qap = E::multi_miller_loop(
            [
                <E::G1Affine as Into<E::G1Prepared>>::into(proof.a),
                prepared_inputs.into_affine().into(),
                proof.c.into(),
            ],
            [
                proof.b.into(),
                pvk.gamma_g2_neg_pc.clone(),
                pvk.delta_g2_neg_pc.clone(),
            ],
        );

        let test = E::final_exponentiation(qap).ok_or(SynthesisError::UnexpectedIdentity)?;

        Ok(test.0 == pvk.alpha_g1_beta_g2)
    }
}
