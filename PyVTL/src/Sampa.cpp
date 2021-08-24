// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2020, Peter Birkholz, Dresden, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
// ****************************************************************************

#include "Sampa.h"

const string Sampa::PHONEME[NUM_PHONEMES] =
{
  // Vowels
  "a:", "e:", "i:", "o:", "u:", "E:", "2:", "y:",
  "a", "e", "i", "o", "u", "E", "2", "y",
  "I", "O", "U", "9", "Y", "@", "6",
  // Diphthongs
  "aI", "aU", "OY",
  "i:6", "i6", "I6", "y:6", "y6", "Y6", "e:6", "e6", "E6", "E:6", 
  "2:6", "26", "96", "a:6", "a6", "u:6", "u6", "U6", "o:6", "o6", "O6",
  // Plosives
  "?", "p", "b", "t", "d", "k", "g",
  // Fricatives
  "f", "v", "T", "D", "s", "z", "S", "Z", "C", "j", "x", "r", "R", "h",
  // Affricates
  "pf", "ts", "tS", "dZ",
  // Nasals and lateral
  "m", "n", "N", "l"
};

// ****************************************************************************
/// Returns the index for a given SAMPA symbol.
// ****************************************************************************

int Sampa::getIndex(const string &symbol)
{
  int index = -1;
  int i;

  for (i=0; (i < NUM_PHONEMES) && (index == -1); i++)
  {
    if (symbol == PHONEME[i]) { index = i; }    
  }

  return index;
}


// ****************************************************************************
/// Is symbol a valid SAMPA phoneme ?
// ****************************************************************************

bool Sampa::isPhoneme(const string &symbol)
{
  return (getIndex(symbol) != -1);
}


// ****************************************************************************
// ****************************************************************************

bool Sampa::isVowel(const string &symbol)
{
  int index = getIndex(symbol);
  if ((index >= 0) && (index < NUM_VOWELS)) 
  { 
    return true; 
  } 
  else 
  { 
    return false; 
  }
}


// ****************************************************************************
// ****************************************************************************

bool Sampa::isDiphthong(const string& symbol)
{
  int index = getIndex(symbol);
  if ((index >= NUM_VOWELS) && (index < NUM_VOWELS + NUM_DIPHTHONGS)) 
  { 
    return true; 
  }
  else 
  { 
    return false; 
  }
}
// ****************************************************************************
// ****************************************************************************

bool Sampa::isConsonant(const string& symbol)
{
  int index = getIndex(symbol);
  if ((index >= NUM_VOWELS + NUM_DIPHTHONGS) && (index < NUM_PHONEMES)) 
  { 
    return true; 
  } 
  else 
  { 
    return false; 
  }
}


// ****************************************************************************
// ****************************************************************************

bool Sampa::isFricative(const string &symbol)
{
  if ((symbol == "f") || (symbol == "v") || 
      (symbol == "T") || (symbol == "D") ||
      (symbol == "s") || (symbol == "z") ||
      (symbol == "S") || (symbol == "Z") ||
      (symbol == "C") || (symbol == "j") || 
      (symbol == "x") || (symbol == "R") || 
      (symbol == "r") || (symbol == "h")) 
  { 
    return true; 
  }

  return false;
}


// ****************************************************************************
// ****************************************************************************

bool Sampa::isPlosive(const string &symbol)
{
  if ((symbol == "b") || (symbol == "p") || 
      (symbol == "d") || (symbol == "t") || 
      (symbol == "g") || (symbol == "k") ||
      (symbol == "?"))
  { 
    return true; 
  }

  return false;
}


// ****************************************************************************
// ****************************************************************************

bool Sampa::isNasal(const string &symbol)
{
  if ((symbol == "m") || (symbol == "n") || (symbol == "N")) 
  { 
    return true; 
  }

  return false;
}


// ****************************************************************************
// ****************************************************************************

bool Sampa::isLateral(const string &symbol)
{
  if (symbol == "l") 
  { 
    return true; 
  }

  return false;
}


// ****************************************************************************
// ****************************************************************************

bool Sampa::isAlveolar(const string& symbol)
{
  if ((symbol == "d") || (symbol == "t") || 
      (symbol == "n") || (symbol == "l") ||
      (symbol == "s") || (symbol == "z"))
  {
    return true;
  }

  return false;
}


// ****************************************************************************
// ****************************************************************************

bool Sampa::isVoiced(const string& symbol)
{
  int index = getIndex(symbol);
  if ((index >= 0) && (index < NUM_VOWELS + NUM_DIPHTHONGS)) 
  { 
    return true; 
  }

  if ((symbol == "b") || (symbol == "d") || (symbol == "g") ||
    (symbol == "v") || (symbol == "D") || (symbol == "z") || 
    (symbol == "Z") || (symbol == "j") || (symbol == "R") || (symbol == "r") ||
    (symbol == "m") || (symbol == "n") || (symbol == "N") || 
    (symbol == "l")) 
  {
    return true;
  }

  return false;
}

// ****************************************************************************
