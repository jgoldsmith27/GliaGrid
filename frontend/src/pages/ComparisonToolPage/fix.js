const fs = require('fs');
const path = require('path');

const filePath = path.join(__dirname, 'ComparisonToolPage.tsx');
let content = fs.readFileSync(filePath, 'utf8');

// Replace the end of the file to add the missing closing brace
const fixed = content.replace(/(\s*<\/Paper>\s*\);)\s*(\n\s*export default ComparisonToolPage;)/s, '$1\n}; // Fixed missing closing brace$2');

fs.writeFileSync(filePath, fixed);
console.log('Fixed the file'); 