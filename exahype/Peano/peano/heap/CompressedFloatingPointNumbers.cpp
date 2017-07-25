#include "peano/heap/CompressedFloatingPointNumbers.h"
#include "tarch/Assertions.h"
#include "tarch/compiler/CompilerSpecificSettings.h"


#include <bitset>



int peano::heap::findMostAgressiveCompression(
  double        values[],
  int           count,
  double        maxError
) {
  assertion(count>0);
  int result = 0;
  for (int i=0; i<count; i++) {
    result = std::max(result, findMostAgressiveCompression(values[i],maxError));
  }
  assertion(result>0);
  return result;
}


#ifdef CompilerICC
int peano::heap::findMostAgressiveCompression(
#else
int __attribute__((optimize("O0"))) peano::heap::findMostAgressiveCompression(
#endif
  double        value,
  double        maxError
) {
  assertion(value==value);
  assertion1( sizeof(long int)*8>=64, sizeof(long int) );

  // We may not use 7 though we use seven out of eight bits for the first byte.
  // If we shift by seven, we can end up with the highest byte set for
  // 0.0155759, e.g.
  int shiftExponent = 6;

  const long int sign = value < 0.0 ? -1 : 1;

  if (sign<0) {
    value = -value;
  }

  int           integerExponent;
  const double  significand          = std::frexp(value , &integerExponent);

  int    usedBytes = 0;
  double error     = std::numeric_limits<double>::max();
  while (
    usedBytes<6
    &&
    error > maxError
  ) {
    usedBytes++;

    const double shiftMantissa    = std::pow( 2.0,shiftExponent );

    double exponent  = static_cast<char>( integerExponent-shiftExponent );
    long int mantissa  = static_cast<long int>( std::round(significand*shiftMantissa) );

    error     = std::abs( std::ldexp(mantissa,exponent) - value );

    shiftExponent+=8;
  }

  assertion(usedBytes>=1);
  assertion(usedBytes<=7);

  return usedBytes;
}



#ifdef CompilerICC
void peano::heap::decomposeIntoEightVariants(
#else
void __attribute__((optimize("O0"))) peano::heap::decomposeIntoEightVariants(
#endif
  double        value,
  char          exponent[8],
  long int      mantissa[8],
  double        error[8]
) {
  assertion(value==value);
  assertion1( sizeof(long int)*8>=64, sizeof(long int) );

  // We may not use 7 though we use seven out of eight bits for the first byte.
  // If we shift by seven, we can end up with the highest byte set for
  // 0.0155759, e.g.
  int shiftExponent = 6;

  const long int sign = value < 0.0 ? -1 : 1;

  if (sign<0) {
    value = -value;
  }

  int           integerExponent;
  const double  significand          = std::frexp(value , &integerExponent);

  for (int i=0; i<8; i++) {
    const double shiftMantissa    = std::pow( 2.0,shiftExponent );

    exponent[i]  = static_cast<char>( integerExponent-shiftExponent );
    mantissa[i]  = static_cast<long int>( std::round(significand*shiftMantissa) );
    error[i]     = std::abs( std::ldexp(mantissa[i],exponent[i]) - value );

    assertion5( mantissa[i]>=0, value, mantissa[i], exponent[i], error[i], sign );

    std::bitset<64>*  mantissaAsBitset = reinterpret_cast<std::bitset<64>*>( &(mantissa[i]) );

    #ifdef Asserts
    for (int j=(i+1)*8-1; j<64; j++) {
      assertion9(
        !(*mantissaAsBitset)[j],
        *mantissaAsBitset, value, static_cast<int>( exponent[i] ), mantissa[i], error[i], i, j, significand, integerExponent
      );
    }
    #endif

    if (sign<0) {
      assertion ( (*mantissaAsBitset)[ (i+1)*8-1 ]==false );
      mantissaAsBitset->flip( (i+1)*8-1 );
    }

    shiftExponent+=8;
  }
}


#ifdef CompilerICC
void peano::heap::decompose(
#else
void __attribute__((optimize("O0"))) peano::heap::decompose(
#endif
  double        value,
  char&         exponent,
  long int&     mantissa,
  int           bytesForMantissa
) {
  assertion(value==value);
  assertion1( sizeof(long int)*8>=64, sizeof(long int) );

  // We may not use 7 though we use seven out of eight bits for the first byte.
  // If we shift by seven, we can end up with the highest byte set for
  // 0.0155759, e.g.
  int shiftExponent = 6;

  const long int sign = value < 0.0 ? -1 : 1;

  if (sign<0) {
    value = -value;
  }

  int           integerExponent;
  const double  significand          = std::frexp(value , &integerExponent);

  for (int i=0; i<bytesForMantissa-1; i++) {
    shiftExponent+=8;
  }

  const double shiftMantissa    = std::pow( 2.0,shiftExponent );

  exponent  = static_cast<char>( integerExponent-shiftExponent );
  mantissa  = static_cast<long int>( std::round(significand*shiftMantissa) );

  assertion4( mantissa>=0, value, mantissa, exponent, sign );

  std::bitset<64>*  mantissaAsBitset = reinterpret_cast<std::bitset<64>*>( &(mantissa) );

  if (sign<0) {
    assertion ( (*mantissaAsBitset)[ (bytesForMantissa)*8-1 ]==false );
    mantissaAsBitset->flip( (bytesForMantissa)*8-1 );
  }
}



void peano::heap::decomposeIntoFourVariants(
  double   value,
  char     exponent[4],
  int      mantissa[4],
  double   error[4]
) {
  assertion(value==value);
  assertion1( sizeof(int)*8>=32, sizeof(int) );

  // We may not use 7 though we use seven out of eight bits for the first byte.
  // If we shift by seven, we can end up with the highest byte set for
  // 0.0155759, e.g.
  int shiftExponent = 6;

  const long int sign = value < 0.0 ? -1 : 1;

  if (sign<0) {
    value = -value;
  }

  int           integerExponent;
  const double  significand          = std::frexp(value , &integerExponent);

  for (int i=0; i<4; i++) {
    const double shiftMantissa    = std::pow( 2.0,shiftExponent );

    exponent[i]  = static_cast<char>( integerExponent-shiftExponent );
    mantissa[i]  = static_cast<int>( std::round(significand*shiftMantissa) );
    error[i]     = std::abs( std::ldexp(mantissa[i],exponent[i]) - value );

    assertion5( mantissa[i]>=0, value, mantissa[i], exponent[i], error[i], sign );

    std::bitset< sizeof(int)*8 >*  mantissaAsBitset = reinterpret_cast<std::bitset< sizeof(int)*8 >*>( &(mantissa[i]) );

    #ifdef Asserts
    for (int j=(i+1)*8-1; j<static_cast<int>( sizeof(int) )*8; j++) {
      assertion9(
        !(*mantissaAsBitset)[j],
        *mantissaAsBitset, value, static_cast<int>( exponent[i] ), mantissa[i], error[i], i, j, significand, integerExponent
      );
    }
    #endif

    if (sign<0) {
      assertion ( (*mantissaAsBitset)[ (i+1)*8-1 ]==false );
      mantissaAsBitset->flip( (i+1)*8-1 );
    }

    shiftExponent+=8;
  }
}


double peano::heap::compose(
  char         exponent,
  long int     mantissa,
  int          bytesUsed
) {
  assertion( bytesUsed>=1 );

  std::bitset<64>*  mantissaAsBitset = reinterpret_cast<std::bitset<64>*>( &mantissa );

  if ( (*mantissaAsBitset)[ bytesUsed*8-1 ] ) {
    mantissaAsBitset->flip( bytesUsed*8-1 );
    mantissa = -mantissa;
  }

  const double doubleMantissa = mantissa;
  const int    intExponent    = exponent;

  return std::ldexp(doubleMantissa,intExponent);
}


