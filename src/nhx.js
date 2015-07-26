/**
 * NHX format parser in JavaScript.
 *
 * Copyright (c) Zorji 2014.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * Example tree (from http://en.wikipedia.org/wiki/Newick_format):
 *
 * +--0.1--A
 * F-----0.2-----B            +-------0.3----C
 * +------------------0.5-----E
 *                            +---------0.4------D
 *
 * Newick format:
 * (A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;
 *
 * Converted to JSON:
 * {
 *   name: "F",
 *   branchset: [
 *     {name: "A", length: 0.1},
 *     {name: "B", length: 0.2},
 *     {
 *       name: "E",
 *       length: 0.5,
 *       branchset: [
 *         {name: "C", length: 0.3},
 *         {name: "D", length: 0.4}
 *       ]
 *     }
 *   ]
 * }
 *
 * ---------------------------------------------------------------------------------------
 * NHX (from http://goby.compbio.cs.cmu.edu/Notung/2.6documentation/Manual-2.6.master014.html)
 *
 * NHX File Format is based on the Newick file format, but embeds additional
 * information about each node in the tree in the comment fields, as follows:
 *
 * [&&NHX:TagID1=value1:TagID2=value2]
 *
 * where TagID1 and TagID2 can specify bootstrap values, species labels, or
 * duplication information. This example has two tags, but NHX comments can
 * have one or more tags. Trees saved in NHX file format include information
 * produced by a reconciliation, including duplications and species labels, but
 * do not record any visual annotations made in Notung. Nor do they record the
 * species tree with which the gene tree was reconciled.
 *
 * NOTE: The NHX format is case-sensitive.
 * More information about NHX format, including a complete list of tags used
 * in comment fields, can be obtained at:
 *
 * http://www.genetics.wustl.edu/eddy/forester/NHX.html.
 *
 * NHX format:
 * (A:0.1[&&NHX:gn=10],B:0.2[&&NHX:gn=10],(C:0.3[&&NHX:gn=10],D:0.4[&&NHX:gn=10])E:0.5[&&NHX:gn=10])F[&&NHX:gn=10];
 *
 * Converted to JSON:
 * {
 *   name: "F",
 *   length: 0,
 *   gn: '10',
 *   branchset: [
 *     {name: "A", length: 0.1, gn: '10'},
 *     {name: "B", length: 0.2, gn: '10'},
 *     {
 *       name: "E",
 *       length: 0.5,
 *       gn: '10',
 *       branchset: [
 *         {name: "C", length: 0.3, gn: '10'},
 *         {name: "D", length: 0.4, gn: '10'}
 *       ]
 *     }
 *   ]
 * }
 *
 */
(function (exports) {
    exports.parse = function (str) {

        function tokenize(str) {
            var current = {
                    branchset: []
                },
                textBegin = 0,
                inNhxAttr = false,
                prevNode = null;

            for (var i = 0; i < str.length; i++) {
                var char = str[i];
                switch (char) {
                    case '(':
                        var parent = current;
                        current = {
                            branchset: []
                        };
                        parent.branchset.push(current);
                        current.parent = parent;
                        textBegin = i + 1;
                        break;
                    case ')':
                        prevNode = current;
                        current = current.parent;
                        textBegin = i + 1;
                        break;
                    case ',':
                        if (inNhxAttr) continue; // ignores comma appear within NHX attributes
                        var text = str.substr(textBegin, i - textBegin);
                        /*
                         * If it is a ')' right before the {text}, then it is
                         * the text of previous node; otherwise it is reading
                         * the text of the sub-nodes of the current node.
                         */
                        var isCharBeforeTextCloseBracket = str.substr(textBegin - 1, 1) == ')';
                        if (isCharBeforeTextCloseBracket) {
                            prevNode.text = text;
                        } else {
                            current.branchset.push({text: text});
                        }

                        textBegin = i + 1;
                        break;
                    case '[':
                        inNhxAttr = true;
                        break;
                    case ']':
                        inNhxAttr = false;
                }
            }

            return current.branchset[0];
        }

        /**
         * Deletes the parent attribute recursively which is only needed in tokenize stage
         *
         * @see tokenize
         * @param {Object} node
         * @returns {Object}
         */
        function deleteParentAttr(node) {
            if (node.hasOwnProperty('branchset')) {
                node.branchset.forEach(deleteParentAttr);
                delete node.parent;
            }
            return node;
        }

        /**
         * input:
         *   A:0.1
         * output:
         *   { name: 'A', length: 0.1 }
         *
         * input:
         *   A
         * output:
         *   { name: 'A', length: 0 }
         *
         * input:
         *  (empty string)
         * output:
         *   { name: '', length: 0 }
         *
         * @param {Object} node
         * @param {string} str
         */
        function setNameAndLength(node, str) {
            if (str.indexOf(':') != -1) {
                var parts = str.split(':');
                node.name = parts[0];
                node.length = parseFloat(parts[1]);
            }
            else if (str.length) {
                node.name = str;
                node.length = 0;
            }
            else {
                node.name = '';
                node.length = 0;
            }
        }

        /**
         * input:
         *   gn=0.2948622:b=6, 7, 10
         * output:
         *   { gn: '0.2948622', b: '6, 7, 10' }
         *
         * @param {Object} node
         * @param {string} text
         */
        function setNhxAttrs(node, text) {
            if (text.length == 0) return;
            foreach(text.split(':'), function (kv) {
                kv = kv.split('=');
                node[kv[0]] = kv[1];
            });
        }

        /**
         * before:
         *   { text: A:0.1[&&NHX:gn=0.2948622:b=6, 7, 10] }
         * after:
         *   { name: 'A', length: 0.1, gn: '0.2948622', b: '6, 7, 10' }
         *
         * @param {Object} node
         * @returns {Object}
         */
        function translateNode(node) {
            if (node.hasOwnProperty('branchset')) {
                foreach(node.branchset, translateNode);
            }
            if (node.hasOwnProperty('text')) {
                var text = node.text.replace(/]$/g, '');
                delete node.text;
                if (text.indexOf('[&&NHX:') != -1) {
                    var parts = text.split('[&&NHX:');
                    setNameAndLength(node, parts[0]);
                    setNhxAttrs(node, parts[1]);
                }
                else {
                    setNameAndLength(node, text);
                }
            }

            return node;
        }

        function foreach(collection, callback) {
            for (var i = 0; i < collection.length; i++) {
                callback(collection[i]);
            }
        }

        /**
         * Appends a trailing comma on every array item to make parsing a bit easier
         *
         * input:
         *   (A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;
         * output:
         *   (A:0.1,B:0.2,(C:0.3,D:0.4,)E:0.5,)F;
         *
         * @param {string} text
         * @returns {string}
         */
        function appendBoundaries(text) {
            text = text.replace(/\)/g, ',)');
            text = text.replace(/;$/, ',');
            return text;
        }

        /**
         * '(A:0.1[&&NHX:gn=10],B:0.2[&&NHX:gn=10],(C:0.3[&&NHX:gn=10],D:0.4[&&NHX:gn=10])E:0.5[&&NHX:gn=10])F[&&NHX:gn=10];'
         */
        var tmp = str;

        /**
         * '(A:0.1[&&NHX:gn=10],B:0.2[&&NHX:gn=10],(C:0.3[&&NHX:gn=10],D:0.4[&&NHX:gn=10],)E:0.5[&&NHX:gn=10],)F[&&NHX:gn=10];'
         */
        tmp = appendBoundaries(tmp);

        /**
         * { branchset:
         *    [ { branchset:
         *         [ { text: 'A:0.1[&&NHX:gn=10]' },
         *           { text: 'B:0.2[&&NHX:gn=10]' },
         *           { branchset:
         *              [ { text: 'C:0.3[&&NHX:gn=10]' },
         *                { text: 'D:0.4[&&NHX:gn=10]' } ],
         *             parent: [Circular],
         *             text: 'E:0.5[&&NHX:gn=10]' } ],
         *        parent: [Circular] } ] }
         */
        tmp = tokenize(tmp);

        /**
         * { branchset:
         *    [ { branchset:
         *         [ { text: 'A:0.1[&&NHX:gn=10]' },
         *           { text: 'B:0.2[&&NHX:gn=10]' },
         *           { branchset:
         *              [ { text: 'C:0.3[&&NHX:gn=10]' },
         *                { text: 'D:0.4[&&NHX:gn=10]' } ],
         *             text: 'E:0.5[&&NHX:gn=10]' } ] } ] }
         */
        tmp = deleteParentAttr(tmp);

        /**
         * { branchset:
         *    [ { branchset:
         *         [ { name: 'A', length: 0.1, gn: '10' },
         *           { name: 'B', length: 0.2, gn: '10' },
         *           { branchset:
         *              [ { name: 'C', length: 0.3, gn: '10' },
         *                { name: 'D', length: 0.4, gn: '10' } ],
         *             name: 'E',
         *             length: 0.5,
         *             gn: '10' } ] } ] }
         */
        tmp = translateNode(tmp);

        return tmp;
    };
})(
    typeof exports !== 'undefined' ? exports : this.NHX = {}
);
