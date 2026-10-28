// The parts of a report page: a titled panel, a part of one, and a block of figures.

function reportPanel(root, title, note = null) {
    const section = document.createElement("section");
    section.className = "panel";
    const heading = document.createElement("h2");
    heading.textContent = title;
    section.appendChild(heading);
    if (note) {
        const text = document.createElement("p");
        text.className = "note";
        text.textContent = note;
        section.appendChild(text);
    }
    root.appendChild(section);
    return section;
}

function panelPart(parent, id, className = null) {
    const part = document.createElement("div");
    if (id) {
        part.id = id;
    }
    if (className) {
        part.className = className;
    }
    parent.appendChild(part);
    return part;
}

// Figures in columns, each column a list of [label, value, id] (id optional, for updating the
// figure later), laid out as given so that related figures sit side by side.
function statGrid(columns) {
    const grid = document.createElement("dl");
    grid.className = "stat-grid";
    grid.style.gridTemplateColumns = `repeat(${columns.length}, minmax(0, 240px))`;
    columns.forEach((entries, column) => {
        entries.forEach(([label, value, id], row) => {
            const item = document.createElement("div");
            item.style.gridColumn = String(column + 1);
            item.style.gridRow = String(row + 1);
            const term = document.createElement("dt");
            term.textContent = label;
            const figure = document.createElement("dd");
            figure.textContent = value;
            if (id) {
                figure.id = id;
            }
            item.append(term, figure);
            grid.appendChild(item);
        });
    });
    return grid;
}
