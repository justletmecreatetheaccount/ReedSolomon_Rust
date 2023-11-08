use std::{str::from_utf8, num::ParseIntError};

struct ReedSolomon {
    pol: u8,
    message: Vec<u8>,
    interpolation_points: Vec<u8>,
}   

impl ReedSolomon {
    
fn byteStr_to_byte(s: &str) -> u8 {
    let mut byte: u8 = 0;
    for c in s.chars() {
        byte = byte << 1;
        if c == '1' {
            byte = byte | 1;
        }
    }
    return byte
}

fn bytes_to_str(byte_vec: Vec<u8>) -> Result<String, std::string::FromUtf8Error> {
    return String::from_utf8(byte_vec);
}

fn multiply_u8(mut a: u8, mut b: u8, pol: u8) -> u8 {
    let mut mult: u8 = 0;
    while b != 0{
        if (b & 1) == 1 {
            mult = mult ^ a;
        }
        if (a & 0b1000_0000) == 0b1000_0000 {
            a = a << 1;
            a = a ^ pol;
        } else {
            a = a << 1;
        }
        b = b >> 1;
    }
    return mult
}

fn evaluate_u8(polynome: &Vec<u8>, a: u8, pol: u8) -> u8 {
    let mut result = polynome.last().cloned().unwrap();
    for i in (0..polynome.len() - 1).rev() {
        result = ReedSolomon::multiply_u8(result, a, pol);
        result = result ^ polynome[i];
    }
    return result;
}

fn invert_u8(a: u8, pol: u8) -> u8 {
    let a2 = ReedSolomon::multiply_u8(a, a, pol);
    let a4 = ReedSolomon::multiply_u8(a2, a2, pol);
    let a8 = ReedSolomon::multiply_u8(a4, a4, pol);
    let a16 = ReedSolomon::multiply_u8(a8, a8, pol);
    let a32 = ReedSolomon::multiply_u8(a16, a16, pol);
    let a64 = ReedSolomon::multiply_u8(a32, a32, pol);
    let a128 = ReedSolomon::multiply_u8(a64, a64, pol);
    return ReedSolomon::multiply_u8(a128, ReedSolomon::multiply_u8(a64, ReedSolomon::multiply_u8(a32, 
        ReedSolomon::multiply_u8(a16, ReedSolomon::multiply_u8(a8, ReedSolomon::multiply_u8(a4, a2, pol), pol), pol), pol), pol), pol);
}

fn gaussian_elimination(y: &Vec<u8>, I: &Vec<u8>, message_len: usize, pol: u8) -> Vec<u8> {
    //create extended matrix
    let mut extended_matrix: Vec<Vec<u8>> = Vec::new();
    for j in 0..message_len {
        extended_matrix.push(Vec::new());
        extended_matrix[j].push(1);
        for i in 1..message_len {
            let prev_val = extended_matrix[j][i-1];
            extended_matrix[j].push(ReedSolomon::multiply_u8(y[j], prev_val, pol));
        }
        extended_matrix[j].push(I[j]);
    }
    //end create extended matrix

    //gaussian elimination
    for i in 0..message_len {
        for j in 0..message_len {
            if i != j {
                let r = ReedSolomon::multiply_u8(extended_matrix[j][i], ReedSolomon::invert_u8(extended_matrix[i][i], pol), pol);
                for m in 0..message_len + 1 {
                    extended_matrix[j][m] = extended_matrix[j][m] ^ ReedSolomon::multiply_u8(r, extended_matrix[i][m], pol);
                }
            }
        }
    }

    let mut result: Vec<u8> = Vec::new();

    for i in 0..message_len {
        result.push(ReedSolomon::multiply_u8(extended_matrix[i][message_len], ReedSolomon::invert_u8(extended_matrix[i][i], pol), pol));
    }

    return result;
}

fn reed_solomon_encode(message: &Vec<u8>, interpolation_points: &Vec<u8>, pol: u8) -> Vec<u8> {
    let mut result: Vec<u8> = Vec::new();
    for x in interpolation_points {
        result.push(ReedSolomon::evaluate_u8(&message, *x, pol));
    }
    return result;
}

fn reed_solomon_decode(byte_vec: &Vec<u8>, message_len: usize, interpolation_points: &Vec<u8>, pol: u8, corrupted_bytes_pos: &Vec<u8>) -> Result<Vec<u8>, bool> {
    let mut adjusted_byte_vec: Vec<u8> = Vec::new();
    let mut adjusted_interpolation_points: Vec<u8> = Vec::new();
    for i in 0..byte_vec.len() {
        if !corrupted_bytes_pos.contains(&(i as u8)) {
            adjusted_byte_vec.push(byte_vec[i]);
            adjusted_interpolation_points.push(interpolation_points[i]);
        }
    }

    if adjusted_byte_vec.len() < message_len {
        return Err(false);
    }

    return Ok(ReedSolomon::gaussian_elimination(&adjusted_interpolation_points, &adjusted_byte_vec, message_len, pol));
}

}
fn main() {
    let test = "TEST String for Encoding and decoding";
    let mon_reed = ReedSolomon{
        pol: 0b01001101,
        message: test.as_bytes().to_vec(),
        interpolation_points: (0..test.len()+10).map(|x| x as u8).collect(),
        }; 

    let encoded_message = ReedSolomon::reed_solomon_encode(&mon_reed.message, &mon_reed.interpolation_points, mon_reed.pol);
    let test_result = ReedSolomon::reed_solomon_decode(
        &encoded_message,
        test.len(),
        &mon_reed.interpolation_points,
        mon_reed.pol,
        &vec![3, 7, 8, 9, 10]);

    println!("{:?}", ReedSolomon::bytes_to_str(test_result.unwrap()));
}
