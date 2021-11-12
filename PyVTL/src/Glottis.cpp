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

#include "Glottis.h"

#include <iomanip>
#include <cstdio>

// This constant should have the same value as the same constant in Tube.h/cpp
const double Glottis::DEFAULT_ASPIRATION_STRENGTH_DB = -40.0;

// ****************************************************************************
/// Default implementation.
// ****************************************************************************

double Glottis::getAspirationStrength_dB()
{
  // When this function was not overwritten by a derived class,
  // just return the default value.
  return DEFAULT_ASPIRATION_STRENGTH_DB;
}


// ****************************************************************************
/// Returns the shape with the given name, or NULL, if a shape with that name
/// is not in the list.
// ****************************************************************************

Glottis::Shape *Glottis::getShape(const string &name)
{
  int i;
  Shape *s = NULL;
  int numShapes = (int)shape.size();

  for (i=0; (i < numShapes) && (s == NULL); i++)
  {
    if (shape[i].name == name)
    {
      s = &shape[i];
    }
  }

  return s;
}


// ****************************************************************************
// Returns the index of the glottis shape with the given name.
// ****************************************************************************

int Glottis::getShapeIndex(const string& name)
{
    int k, index = -1;
    int numShapes = (int)shape.size();
    for (k = 0; (k < numShapes) && (index == -1); k++)
    {
        if (shape[k].name == name)
        {
            index = k;
        }
    }
    return index;
}

// ****************************************************************************
/// Have any changes to the shapes been made since the last saving?
// ****************************************************************************

bool Glottis::hasUnsavedChanges()
{
  int i, k;
  bool changes = false;

  if ((savedState.staticParam.size() != staticParam.size()) ||
      (savedState.shape.size() != shape.size()))
  {
    return true;
  }

  // Any changes in the static parameters ?

  for (i=0; i < (int)staticParam.size(); i++)
  {
    if (staticParam[i].x != savedState.staticParam[i])
    {
      return true;
    }
  }

  // Any changes in the shapes ?

  for (i=0; i < (int)shape.size(); i++)
  {
    if ((savedState.shape[i].controlParam.size() != controlParam.size()) ||
        (savedState.shape[i].name != shape[i].name))
    {
      return true;
    }

    for (k=0; k < (int)controlParam.size(); k++)
    {
      if (savedState.shape[i].controlParam[k] != shape[i].controlParam[k])
      {
        return true;
      }
    }
  }

  // Return "no changes".
  return false;
}


// ****************************************************************************
/// Declare the current settings as the saved state.
// ****************************************************************************

void Glottis::clearUnsavedChanges()
{
  int i;

  savedState.shape = shape;

  savedState.staticParam.resize( staticParam.size() );
  for (i=0; i < (int)staticParam.size(); i++)
  {
    savedState.staticParam[i] = staticParam[i].x;
  }
}


// ****************************************************************************
/// Writes the glottis data to an XML string.
// ****************************************************************************

bool Glottis::writeToXml(ostream &os, int initialIndent, bool isSelected)
{
  int i, k;
  char st[1024];
  int indent = initialIndent;

  // ****************************************************************
  // Open the glottis element.
  // ****************************************************************

  os << string(indent, ' ') << "<glottis_model type=\"" << getName() << "\" selected=\"" << isSelected << "\">" << endl;
  indent+= 2;

  // ****************************************************************
  // Write the static parameters.
  // ****************************************************************

  os << string(indent, ' ') << "<static_params>" << endl;
  indent+= 2;

  for (i=0; i < (int)staticParam.size(); i++)
  {
    sprintf(st, "<param index=\"%d\" name=\"%s\" description=\"%s\" unit=\"%s\" min=\"%f\" max=\"%f\" neutral=\"%f\"/>",
      i,
      staticParam[i].name.c_str(),
      staticParam[i].description.c_str(),
      staticParam[i].cgsUnit.c_str(),
      staticParam[i].min,
      staticParam[i].max,
      staticParam[i].neutral);

    os << string(indent, ' ') << st << endl;
  }

  indent-= 2;
  os << string(indent, ' ') << "</static_params>" << endl;

  // ****************************************************************
  // Write the control parameters.
  // ****************************************************************

  os << string(indent, ' ') << "<control_params>" << endl;
  indent+= 2;

  for (i=0; i < (int)controlParam.size(); i++)
  {
    sprintf(st, "<param index=\"%d\" name=\"%s\" description=\"%s\" unit=\"%s\" min=\"%f\" max=\"%f\" neutral=\"%f\"/>",
      i,
      controlParam[i].name.c_str(),
      controlParam[i].description.c_str(),
      controlParam[i].cgsUnit.c_str(),
      controlParam[i].min,
      controlParam[i].max,
      controlParam[i].neutral);

    os << string(indent, ' ') << st << endl;
  }

  indent-= 2;
  os << string(indent, ' ') << "</control_params>" << endl;

  // ****************************************************************
  // Write the shapes.
  // ****************************************************************

  os << string(indent, ' ') << "<shapes>" << endl;
  indent+= 2;

  for (k=0; k < (int)shape.size(); k++)
  {
    os << string(indent, ' ') << "<shape name=\"" << shape[k].name << "\">" << endl;
    indent+= 2;

    for (i=0; i < (int)controlParam.size(); i++)
    {
      sprintf(st, "<control_param index=\"%d\" value=\"%f\"/>", i, shape[k].controlParam[i]);
      os << string(indent, ' ') << st << endl;
    }

    indent-= 2;
    os << string(indent, ' ') << "</shape>" << endl;
  }

  indent-= 2;
  os << string(indent, ' ') << "</shapes>" << endl;

  // ****************************************************************
  // Close the glottis element.
  // ****************************************************************

  indent-= 2;
  os << string(indent, ' ') << "</glottis_model>" << endl;

  // Declare the current state as saved.
  clearUnsavedChanges();

  return true;
}


// ****************************************************************************
/// Reads the glottis data and shapes from the given XML-structure.
// ****************************************************************************

bool Glottis::readFromXml(XmlNode &rootNode)
{
  int i, k, m, n;
  XmlNode *node;
  
  string name;
  double value;
  double min;
  double max;
  double neutral;

  // ****************************************************************
  // Read the static parameter values.
  // ****************************************************************

  XmlNode *staticParamNode = rootNode.getChildElement("static_params");
  
  if (staticParamNode != NULL)
  {
    for (i=0; i < (int)staticParamNode->childElement.size(); i++)
    {
      node = staticParamNode->childElement[i];
      try
      {
          name = node->getAttributeString("name");
          min = node->getAttributeDouble("min");
          max = node->getAttributeDouble("max");
          neutral = node->getAttributeDouble("neutral");
      }
      catch (const std::exception & e)
      {
          printf("Error: %s\n", e.what());
          return false;
      }
      n = -1;
      for (m = 0; (m < staticParam.size()) && (n == -1); m++)
      {
          if (name == staticParam[m].name)
          {
              n = m;
          }
      }
      if (n != -1)
      {
          staticParam[n].min = min;
          staticParam[n].max = max;
          staticParam[n].neutral = neutral;
          staticParam[n].x = neutral;
      }
      else
      {
          printf("Error: Could not read definition of static parameters of selected glottis model! One of the parameter names must be incorrect.\n");
          return false;
      }
    }
  }

  // ****************************************************************
  // Read the control parameter values.
  // ****************************************************************

  XmlNode *controlParamNode = rootNode.getChildElement("control_params");

  if (controlParamNode != NULL)
  {
    for (i=0; i < (int)controlParamNode->childElement.size(); i++)
    {
      node = controlParamNode->childElement[i];
      try
      {
          name = node->getAttributeString("name");
          min = node->getAttributeDouble("min");
          max = node->getAttributeDouble("max");
          neutral = node->getAttributeDouble("neutral");
      }
      catch (const std::exception & e)
      {
          printf("Error: %s\n", e.what());
          return false;
      }
      n = -1;
      for (m = 0; (m < controlParam.size()) && (n == -1); m++)
      {
          if (name == controlParam[m].name)
          {
              n = m;
          }
      }
      if (n != -1)
      {
          controlParam[n].min = min;
          controlParam[n].max = max;
          controlParam[n].neutral = neutral;
          controlParam[n].x = neutral;
      }
      else
      {
          printf("Error: Could not read definition of control parameters of selected glottis model! One of the parameter names must be incorrect.\n");
          return false;
      }
    }
  }

  // ****************************************************************
  // Read the shapes.
  // ****************************************************************

  XmlNode *shapesNode = rootNode.getChildElement("shapes");

  if (shapesNode != NULL)
  {
    // Clear all current shapes
    shape.clear();

    // Run through all shapes
    int numshapes = shapesNode->numChildElements("shape");
    for (i=0; i < numshapes; i++)
    {
      XmlNode *shapeNode = shapesNode->getChildElement("shape", i);
      Shape s;
      s.name = shapeNode->getAttributeString("name");

      // Init default values for the shape to read in.
      s.controlParam.resize( controlParam.size() );
      for (k=0; k < (int)controlParam.size(); k++)
      {
        s.controlParam[k] = controlParam[k].neutral;
      }

      if (s.name.empty())
      {
        printf("Error: Invalid glottis shape element. The shape name is missing!\n");
        return false;
      }

      // Run through all control parameter values
      int numParams = shapeNode->numChildElements("control_param");

      for (k=0; k < numParams; k++)
      {
        node = shapeNode->getChildElement("control_param", k);
        try
        {
            name = node->getAttributeString("name");
            value = node->getAttributeDouble("value");
        }
        catch (const std::exception& e)
        {
            printf("Error: %s\n", e.what());
            return false;
        }
        n = -1;
        for (m = 0; (m < controlParam.size()) && (n == -1); m++)
        {
            if (name == controlParam[m].name)
            {
                n = m;
            }
        }
        if (n != -1)
        {
            s.controlParam[n] = value;
        }
        else
        {
            printf("Error: Could not read control paramter values in glottis shape: %s !\n", s.name.c_str() );
            return false;
        }
      }

      // Add the read shape to the shapes list.
      shape.push_back(s);
    }
  }

  // ****************************************************************
  // Calculate the geometry with the loaded data.
  // ****************************************************************

  resetMotion();
  calcGeometry();

  // Declare the current state as saved.
  clearUnsavedChanges();

  return true;
}

// ****************************************************************************
/// Print the names and units of the control and derived parameters in a row.
// ****************************************************************************

void Glottis::printParamNames(ostream &os)
{
  int i;

  for (i=0; i < (int)controlParam.size(); i++)
  {
    os << controlParam[i].name << "[" << controlParam[i].cgsUnit << "] ";
  }

  for (i=0; i < (int)derivedParam.size(); i++)
  {
    os << derivedParam[i].name << "[" << derivedParam[i].cgsUnit << "] ";
  }

  os << "glottal_flow[cm^3/s] ";
  os << "P_subglottal[dPa] ";
  os << "P_intraglottal1[dPa] ";
  os << "P_intraglottal2[dPa] ";
  os << "P_supraglottal[dPa] ";
  os << "mouth_flow[cm^3/s] ";
  os << "nostril_flow[cm^3/s] ";
  os << "skin_flow[cm^3/s] ";
  os << "radiated_pressure[dPa] ";
  os << endl;
}


// ****************************************************************************
/// Print the values of the control and derived parameters in a row.
// ****************************************************************************

void Glottis::printParamValues(ostream &os, double glottalFlow_cm3_s, 
  double glottalPressure_dPa[], double mouthFlow_cm3_s, double nostrilFlow_cm3_s,
  double skinFlow_cm3_s, double radiatedPressure_dPa)
{
  int i;

  // Specify how many digits to display after the decimal point.
  os << fixed << setprecision(8);
  
  for (i=0; i < (int)controlParam.size(); i++)
  {
    os << controlParam[i].x << " ";
  }

  for (i=0; i < (int)derivedParam.size(); i++)
  {
    os << derivedParam[i].x << " ";
  }

  // Add the current flow value and the four pressure values as last 
  // items in the row.
  os << glottalFlow_cm3_s << " ";
  os << glottalPressure_dPa[0] << " ";
  os << glottalPressure_dPa[1] << " ";
  os << glottalPressure_dPa[2] << " ";
  os << glottalPressure_dPa[3] << " ";
  os << mouthFlow_cm3_s << " ";
  os << nostrilFlow_cm3_s << " ";
  os << skinFlow_cm3_s << " ";
  os << radiatedPressure_dPa;
  
  os << endl;
}

// ****************************************************************************
/// Restricts the values of all parameters in the given vector.
// ****************************************************************************

void Glottis::restrictParams(vector<Parameter> &p)
{
  int i;

  for (i=0; i < (int)p.size(); i++)
  {
    if (p[i].x < p[i].min)
    {
      p[i].x = p[i].min;
    }

    if (p[i].x > p[i].max)
    {
      p[i].x = p[i].max;
    }
  }
}


// ****************************************************************************
/// Temporarily store (cache) the control parameter values, so that they can 
/// restored later.
// ****************************************************************************

void Glottis::storeControlParams()
{
  hasStoredControlParams = true;
  int i;
  for (i = 0; i < (int)controlParam.size(); i++)
  {
    storedControlParams[i] = controlParam[i].x;
  }
}


// ****************************************************************************
/// Restore the temporarily stored (cached) control parameter values.
// ****************************************************************************

void Glottis::restoreControlParams()
{
  if (hasStoredControlParams == false)
  {
    return;
  }

  int i;
  for (i = 0; i < (int)controlParam.size(); i++)
  {
    controlParam[i].x = storedControlParams[i];
  }
  hasStoredControlParams = false;

  calcGeometry();
}


// ****************************************************************************
