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

#ifndef __XML_NODE_H__
#define __XML_NODE_H__

#include "XmlExceptions.h"

#include <string>
#include <vector>


using namespace std;

// ****************************************************************************
// This file implements a simple XML parser using the STL.
// ****************************************************************************

struct XmlError;
struct XmlAttribute;
class XmlNode;

// ****************************************************************************
// Some public functions.
// Use xmlParseString(...) or xmlParseFile(...) to get the root node of the
// XML-document. After processing, the caller must delete the returned node
// with "delete" to free the memory of the XML-tree again.
// If the return value of the parsing functions is NULL, use xmlPrintErrors(...)
// to print the parsing errors.
// ****************************************************************************

XmlNode *xmlParseString(const string &input, const string &tag, vector<XmlError> *errors = NULL);
XmlNode *xmlParseFile(const string &fileName, const string &tag, vector<XmlError> *errors = NULL);
void xmlPrintErrors(vector<XmlError> &errors);
void xmlTest();

// ****************************************************************************
// Structure for an error detected during parsing.
// ****************************************************************************

struct XmlError
{
  int line;
  int column;
  string text;
};

// ****************************************************************************
// Structure for an attribute of an XML element.
// ****************************************************************************

struct XmlAttribute
{
  string name;
  string value;
};

// ****************************************************************************
/// An XML node.
// ****************************************************************************

class XmlNode
{
  // ****************************************************************
  // Public data.
  // ****************************************************************

public:
  enum NodeType
  {
    ELEMENT,
    TEXT,
    OTHER       ///< declaration-nodes, comments, CDATA-nodes, etc.
  };

  // Common members in all types of nodes.
  XmlNode *parent;
  NodeType type;

  // Exclusively used in ELEMENT nodes.
  string name;
  vector<XmlNode*> child;           ///< All child nodes.
  vector<XmlNode*> childElement;    ///< Only element child nodes.
  vector<XmlAttribute> attribute;

  // Exclusively used in TEXT, COMMENT and OTHER nodes.
  // These are leaf nodes without children or attributes.
  string text;

  // ****************************************************************
  // Public functions.
  // ****************************************************************

public:
  XmlNode(NodeType t, XmlNode *parent);
  ~XmlNode();

  int numChildElements(const string &name);
  XmlNode *getChildElement(const string &name, int index = 0);
  bool hasAttribute(const string &name);
  int getAttributeInt(const string &name);
  double getAttributeDouble(const string &name);
  string getAttributeString(const string &name);

  string toXmlString();

  // ****************************************************************
  // Private functions.
  // ****************************************************************

private:
  void toXmlString(ostream &os, int indent);
};


#endif
