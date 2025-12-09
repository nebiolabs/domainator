const state = {
  files: [],
  tools: [],
  jobs: {},
  lastJobUpdate: 0,
  activeToolId: null,
  collapsedPanels: new Set(),
  categoryConfig: {
    order: [],
    colors: {},
    orderLookup: new Map(),
    colorLookup: new Map(),
    defaultColor: "#f9fafb",
  },
};

const uploadState = {
  files: [],
  suppressInputChange: false,
};

const DOCUMENTATION_URL =
  "https://github.com/nebiolabs/domainator/tree/main/docs/server/README.md";

const PREVIEWABLE_TYPES = new Set([
  "html",
  "htm",
  "svg",
  "txt",
  "text",
  "json",
  "csv",
  "tsv",
  "png",
  "jpg",
  "jpeg",
  "gif",
  "pdf",
]);

function supportsDataTransferConstructor() {
  return typeof DataTransfer === "function";
}

function getUploadFileKey(file) {
  if (!file) {
    return "";
  }
  return [file.name, file.size, file.lastModified].join("::");
}

function syncFileInput() {
  const input = document.getElementById("file-input");
  if (!input) {
    return;
  }

  uploadState.suppressInputChange = true;
  try {
    if (uploadState.files.length && supportsDataTransferConstructor()) {
      const dataTransfer = new DataTransfer();
      uploadState.files.forEach(file => {
        try {
          dataTransfer.items.add(file);
        } catch (error) {
          console.warn("Unable to add file to drop selection.", error);
        }
      });
      input.files = dataTransfer.files;
    } else if (!uploadState.files.length) {
      input.value = "";
    }
  } finally {
    uploadState.suppressInputChange = false;
  }
}

function updateDropZoneSummary() {
  const detail = document.getElementById("drop-zone-detail");
  if (!detail) {
    return;
  }
  if (!uploadState.files.length) {
    detail.textContent = "No files selected yet.";
    return;
  }
  if (uploadState.files.length === 1) {
    detail.textContent = uploadState.files[0].name;
    return;
  }
  const maxPreview = 3;
  const preview = uploadState.files.slice(0, maxPreview).map(file => file.name);
  const remaining = uploadState.files.length - preview.length;
  detail.textContent =
    remaining > 0 ? `${preview.join(", ")} + ${remaining} more` : preview.join(", ");
}

function addFilesToUpload(fileList) {
  if (!fileList || !fileList.length) {
    return 0;
  }
  const incoming = Array.from(fileList);
  if (!incoming.length) {
    return 0;
  }
  const existingKeys = new Set(uploadState.files.map(getUploadFileKey));
  const added = [];
  incoming.forEach(file => {
    const key = getUploadFileKey(file);
    if (key && !existingKeys.has(key)) {
      uploadState.files.push(file);
      existingKeys.add(key);
      added.push(file);
    }
  });
  if (added.length) {
    syncFileInput();
    updateDropZoneSummary();
  }
  return added.length;
}

function clearUploadSelection() {
  uploadState.files = [];
  syncFileInput();
  updateDropZoneSummary();
}

function initializeDropZone() {
  const dropZone = document.getElementById("drop-zone");
  const fileInput = document.getElementById("file-input");
  if (!dropZone || !fileInput) {
    return;
  }

  const preventDefaults = event => {
    event.preventDefault();
    event.stopPropagation();
  };

  const activate = () => dropZone.classList.add("is-active");
  const deactivate = () => dropZone.classList.remove("is-active");

  ["dragenter", "dragover"].forEach(eventName => {
    dropZone.addEventListener(eventName, event => {
      preventDefaults(event);
      activate();
    });
  });

  ["dragleave", "dragend"].forEach(eventName => {
    dropZone.addEventListener(eventName, event => {
      preventDefaults(event);
      deactivate();
    });
  });

  dropZone.addEventListener("drop", event => {
    preventDefaults(event);
    deactivate();
    const files = event.dataTransfer?.files;
    const added = addFilesToUpload(files);
    if (added) {
      setMessage(`${added} file${added === 1 ? "" : "s"} ready for upload.`, "success", 4000);
    }
  });

  dropZone.addEventListener("click", () => fileInput.click());
  dropZone.addEventListener("keydown", event => {
    if (event.key === "Enter" || event.key === " ") {
      event.preventDefault();
      fileInput.click();
    }
  });

  fileInput.addEventListener("change", () => {
    if (uploadState.suppressInputChange) {
      return;
    }
    const added = addFilesToUpload(fileInput.files);
    if (!added) {
      updateDropZoneSummary();
    }
  });

  updateDropZoneSummary();
}

function normalizeCategory(category) {
  if (!category || typeof category !== "string") {
    return "\uffff";
  }
  return category.toLowerCase();
}

function getCategoryOrderIndex(category) {
  if (!category || typeof category !== "string") {
    return state.categoryConfig.order.length;
  }
  const index = state.categoryConfig.orderLookup.get(category.toLowerCase());
  return index !== undefined ? index : state.categoryConfig.order.length;
}

function compareCategoryOrder(categoryA, categoryB) {
  const indexA = getCategoryOrderIndex(categoryA);
  const indexB = getCategoryOrderIndex(categoryB);
  if (indexA !== indexB) {
    return indexA - indexB;
  }
  const normalizedA = normalizeCategory(categoryA);
  const normalizedB = normalizeCategory(categoryB);
  if (normalizedA !== normalizedB) {
    return normalizedA.localeCompare(normalizedB);
  }
  return 0;
}

function normalizeToolName(tool) {
  return (tool?.display_name || tool?.id || "").toLowerCase();
}

function getToolPriority(tool) {
  if (!tool || tool.priority === undefined || tool.priority === null) {
    return null;
  }
  const value = Number(tool.priority);
  return Number.isFinite(value) ? value : null;
}

function compareTools(a, b) {
  const categoryComparison = compareCategoryOrder(a?.category, b?.category);
  if (categoryComparison !== 0) {
    return categoryComparison;
  }

  const priorityA = getToolPriority(a);
  const priorityB = getToolPriority(b);
  const hasPriorityA = priorityA !== null;
  const hasPriorityB = priorityB !== null;
  if (hasPriorityA && hasPriorityB) {
    if (priorityA !== priorityB) {
      return priorityA - priorityB;
    }
  } else if (hasPriorityA) {
    return -1;
  } else if (hasPriorityB) {
    return 1;
  }

  const nameA = normalizeToolName(a);
  const nameB = normalizeToolName(b);
  if (nameA !== nameB) {
    return nameA.localeCompare(nameB);
  }
  const idA = (a?.id || "").toLowerCase();
  const idB = (b?.id || "").toLowerCase();
  return idA.localeCompare(idB);
}

function applyCategoryConfig(config) {
  const defaultColor =
    typeof config?.default_color === "string" && config.default_color.trim()
      ? config.default_color.trim()
      : state.categoryConfig.defaultColor;

  const seenOrder = new Set();
  const order = Array.isArray(config?.order)
    ? config.order
        .map(item => (typeof item === "string" ? item.trim() : ""))
        .filter(item => {
          if (!item) {
            return false;
          }
          const key = item.toLowerCase();
          if (seenOrder.has(key)) {
            return false;
          }
          seenOrder.add(key);
          return true;
        })
    : [];

  const colorsInput = config?.colors && typeof config.colors === "object" ? config.colors : {};
  const colors = {};
  const colorLookup = new Map();
  Object.entries(colorsInput).forEach(([rawName, rawValue]) => {
    if (typeof rawName !== "string") {
      return;
    }
    const name = rawName.trim();
    const value = typeof rawValue === "string" ? rawValue.trim() : "";
    if (!name || !value) {
      return;
    }
    colors[name] = value;
    colorLookup.set(name.toLowerCase(), value);
  });

  const orderLookup = new Map();
  order.forEach((name, index) => {
    orderLookup.set(name.toLowerCase(), index);
  });

  state.categoryConfig = {
    order,
    colors,
    orderLookup,
    colorLookup,
    defaultColor,
  };
}

async function loadCategoryConfig() {
  const url = "/static/tool-category-config.json";
  try {
    const response = await fetch(url, { cache: "no-cache" });
    if (!response.ok) {
      throw new Error(`${response.status} ${response.statusText}`);
    }
    const config = await response.json();
    applyCategoryConfig(config);
  } catch (error) {
    console.warn("Unable to load tool category configuration; using defaults.", error);
    applyCategoryConfig({});
  }
}

function clampColorValue(value) {
  return Math.max(0, Math.min(255, value));
}

function parseHexColor(color) {
  if (typeof color !== "string") {
    return null;
  }
  const hex = color.trim();
  const match = /^#?([0-9a-f]{3}|[0-9a-f]{6})$/i.exec(hex);
  if (!match) {
    return null;
  }
  let value = match[1];
  if (value.length === 3) {
    value = value
      .split("")
      .map(char => char + char)
      .join("");
  }
  const intValue = parseInt(value, 16);
  return {
    r: (intValue >> 16) & 0xff,
    g: (intValue >> 8) & 0xff,
    b: intValue & 0xff,
  };
}

function parseRgbColor(color) {
  if (typeof color !== "string") {
    return null;
  }
  const match = /^rgba?\(([^)]+)\)$/i.exec(color.trim());
  if (!match) {
    return null;
  }
  const parts = match[1]
    .split(",")
    .map(part => part.trim())
    .filter(Boolean);
  if (parts.length < 3) {
    return null;
  }
  const r = clampColorValue(Number.parseFloat(parts[0]));
  const g = clampColorValue(Number.parseFloat(parts[1]));
  const b = clampColorValue(Number.parseFloat(parts[2]));
  if (Number.isNaN(r) || Number.isNaN(g) || Number.isNaN(b)) {
    return null;
  }
  return { r, g, b };
}

function parseColor(color) {
  return parseHexColor(color) || parseRgbColor(color);
}

function rgbToHex({ r, g, b }) {
  const toHex = value => clampColorValue(Math.round(value)).toString(16).padStart(2, "0");
  return `#${toHex(r)}${toHex(g)}${toHex(b)}`;
}

function darkenColor(color, amount = 0.16) {
  const parsed = parseColor(color);
  if (!parsed) {
    return color;
  }
  const factor = 1 - amount;
  return rgbToHex({
    r: parsed.r * factor,
    g: parsed.g * factor,
    b: parsed.b * factor,
  });
}

function getCategoryColors(category) {
  let fill = state.categoryConfig.defaultColor;
  if (typeof category === "string") {
    const configured = state.categoryConfig.colorLookup.get(category.toLowerCase());
    if (configured) {
      fill = configured;
    }
  }
  const border = darkenColor(fill, 0.18);
  return {
    fill,
    border,
  };
}

async function initializeApp() {
  initializeDropZone();
  bindEventHandlers();
  initializePanelControls();
  await loadCategoryConfig();
  await refreshAll();
  startPolling();
}

function getExtension(value) {
  if (!value) {
    return "";
  }
  const lastDot = value.lastIndexOf(".");
  if (lastDot < 0 || lastDot === value.length - 1) {
    return "";
  }
  return value.slice(lastDot + 1).toLowerCase();
}

function isPreviewableFile(name, type) {
  const normalizedType = (type || "").toLowerCase();
  if (normalizedType && PREVIEWABLE_TYPES.has(normalizedType)) {
    return true;
  }
  const extension = getExtension((name || "").toLowerCase());
  if (PREVIEWABLE_TYPES.has(extension)) {
    return true;
  }
  if (!extension && normalizedType === "") {
    return false;
  }
  if (normalizedType.includes("html")) {
    return true;
  }
  if (normalizedType.includes("svg")) {
    return true;
  }
  if (normalizedType.startsWith("image/")) {
    return true;
  }
  return false;
}

let messageTimer = null;
let pollTimer = null;
let isPolling = false;

document.addEventListener("DOMContentLoaded", () => {
  initializeApp().catch(error => {
    console.error("Failed to initialize the Domainator interface.", error);
    setMessage("Failed to initialize the interface. Please reload the page.", "error");
  });
});

function bindEventHandlers() {
  const uploadForm = document.getElementById("upload-form");
  if (uploadForm) {
    uploadForm.addEventListener("submit", handleUploadSubmit);
  }

  const toolFilter = document.getElementById("tool-filter");
  if (toolFilter) {
    toolFilter.addEventListener("input", () => renderTools());
  }

  const panelControls = document.getElementById("panel-controls");
  if (panelControls) {
    panelControls.addEventListener("click", event => {
      const button = event.target.closest("button[data-panel]");
      if (!button) {
        return;
      }
      const panelId = button.getAttribute("data-panel");
      if (!panelId) {
        return;
      }
      togglePanel(panelId);
    });
  }

  const documentationButton = document.getElementById("documentation-button");
  if (documentationButton) {
    documentationButton.addEventListener("click", () => {
      const win = window.open(DOCUMENTATION_URL, "_blank");
      if (win) {
        win.opener = null;
      } else {
        window.location.href = DOCUMENTATION_URL;
      }
    });
  }

  const toolList = document.getElementById("tool-list");
  if (toolList) {
    toolList.addEventListener("click", event => {
      const card = event.target.closest(".tool-card");
      if (!card) {
        return;
      }
      selectTool(card.dataset.toolId);
    });
  }

  const filesBody = document.getElementById("files-body");
  if (filesBody) {
    filesBody.addEventListener("click", event => {
      const button = event.target.closest("button[data-delete]");
      if (!button) {
        return;
      }
      const fileId = button.getAttribute("data-delete");
      if (!fileId) {
        return;
      }
      deleteFile(fileId);
    });
  }

  const jobsBody = document.getElementById("jobs-body");
  if (jobsBody) {
    jobsBody.addEventListener("click", event => {
      const button = event.target.closest("button[data-cancel]");
      if (!button) {
        return;
      }
      const jobId = button.getAttribute("data-cancel");
      if (!jobId) {
        return;
      }
      cancelJob(jobId);
    });
  }
}

async function refreshAll() {
  await Promise.all([refreshFiles(), refreshTools(), refreshJobs()]);
}

async function refreshFiles() {
  try {
    const response = await fetch("/api/files");
    if (!response.ok) {
      throw new Error("Failed to fetch files");
    }
    const data = await response.json();
    state.files = Array.isArray(data)
      ? data.sort((a, b) => new Date(b.uploaded_at) - new Date(a.uploaded_at))
      : [];
    renderFiles();
    if (state.activeToolId) {
      renderToolDetail(getActiveTool());
    }
  } catch (error) {
    setMessage(error.message, "error");
  }
}

async function refreshTools() {
  try {
    const response = await fetch("/api/tools");
    if (!response.ok) {
      throw new Error("Failed to fetch tools");
    }
    const data = await response.json();
    const tools = Array.isArray(data) ? data.slice() : [];
    tools.sort(compareTools);
    state.tools = tools;
    if (!state.tools.some(tool => tool.id === state.activeToolId)) {
      state.activeToolId = null;
    }
    renderTools();
    renderToolDetail(getActiveTool());
  } catch (error) {
    setMessage(error.message, "error");
  }
}

async function refreshJobs() {
  try {
    const response = await fetch("/api/jobs");
    if (!response.ok) {
      throw new Error("Failed to fetch jobs");
    }
    const jobs = await response.json();
    const list = Array.isArray(jobs) ? jobs : [];
    state.jobs = {};
    state.lastJobUpdate = 0;
    updateJobs(list);
  } catch (error) {
    setMessage(error.message, "error");
  }
}

function renderFiles() {
  const body = document.getElementById("files-body");
  if (!body) {
    return;
  }
  body.innerHTML = "";
  if (!state.files.length) {
    const row = document.createElement("tr");
    const cell = document.createElement("td");
    cell.colSpan = 5;
    cell.textContent = "No files uploaded yet.";
    row.appendChild(cell);
    body.appendChild(row);
    return;
  }
  state.files.forEach(file => {
    const row = document.createElement("tr");

    const nameCell = document.createElement("td");
    nameCell.textContent = file.original_name;
    row.appendChild(nameCell);

    const typeCell = document.createElement("td");
    typeCell.textContent = file.type || "";
    row.appendChild(typeCell);

    const sizeCell = document.createElement("td");
    sizeCell.textContent = formatBytes(file.size || 0);
    row.appendChild(sizeCell);

    const uploadedCell = document.createElement("td");
    uploadedCell.textContent = formatIsoTimestamp(file.uploaded_at);
    row.appendChild(uploadedCell);

    const actionCell = document.createElement("td");
    const isPreviewable = isPreviewableFile(file.original_name, file.type);
    if (isPreviewable) {
      const viewLink = document.createElement("a");
      viewLink.href = `/api/files/${file.file_id}/view`;
      viewLink.textContent = "View";
      viewLink.className = "view-link";
      viewLink.target = "_blank";
      viewLink.rel = "noopener";
      actionCell.appendChild(viewLink);
      actionCell.appendChild(document.createTextNode(" "));
    }
    const downloadLink = document.createElement("a");
    downloadLink.href = `/api/files/${file.file_id}/download`;
    downloadLink.textContent = "Download";
    downloadLink.className = "download-link";
    downloadLink.target = "_blank";
    actionCell.appendChild(downloadLink);
    actionCell.appendChild(document.createTextNode(" "));

    const deleteButton = document.createElement("button");
    deleteButton.type = "button";
    deleteButton.textContent = "Delete";
    deleteButton.setAttribute("data-delete", file.file_id);
    deleteButton.className = "delete-button";
    actionCell.appendChild(deleteButton);
    row.appendChild(actionCell);

    body.appendChild(row);
  });
}

function renderTools() {
  const container = document.getElementById("tool-list");
  if (!container) {
    return;
  }
  container.innerHTML = "";
  const filterValue = (document.getElementById("tool-filter")?.value || "").trim().toLowerCase();

  const filtered = state.tools.filter(tool => {
    if (!filterValue) {
      return true;
    }
    const haystack = [tool.id, tool.display_name, tool.category]
      .filter(Boolean)
      .join(" ")
      .toLowerCase();
    return haystack.includes(filterValue);
  });

  if (!filtered.length) {
    const empty = document.createElement("p");
    empty.textContent = state.tools.length ? "No tools matched that filter." : "No tools have been registered.";
    container.appendChild(empty);
    return;
  }

  filtered.forEach(tool => {
    const card = document.createElement("div");
    card.className = "tool-card";
    card.dataset.toolId = tool.id;
    if (tool.id === state.activeToolId) {
      card.classList.add("active");
    }

    const { fill, border } = getCategoryColors(tool.category);
    card.style.setProperty("--tool-card-fill", fill);
    card.style.setProperty("--tool-card-border", border);

    const title = document.createElement("h3");
    title.textContent = tool.display_name || tool.id;
    card.appendChild(title);

    if (tool.category) {
      const category = document.createElement("div");
      category.className = "category";
      category.textContent = tool.category;
      card.appendChild(category);
    }

    if (tool.description) {
      const desc = document.createElement("p");
      const firstLine = tool.description.split(/\r?\n/, 1)[0].trim();
      if (firstLine) {
        desc.textContent = firstLine;
        card.appendChild(desc);
      }
    }

    container.appendChild(card);
  });
}

function selectTool(toolId) {
  if (!toolId) {
    state.activeToolId = null;
    renderTools();
    renderToolDetail(null);
    return;
  }
  openPanel("tool-detail-panel");
  scrollWindowToTop();
  if (state.activeToolId === toolId) {
    return;
  }
  state.activeToolId = toolId;
  renderTools();
  renderToolDetail(getActiveTool());
}

function getActiveTool() {
  if (!state.activeToolId) {
    return null;
  }
  return state.tools.find(tool => tool.id === state.activeToolId) || null;
}

function renderToolDetail(tool) {
  const container = document.getElementById("tool-detail");
  if (!container) {
    return;
  }
  container.innerHTML = "";

  if (!tool) {
    const message = document.createElement("p");
    message.textContent = "Select a tool to see parameters.";
    container.appendChild(message);
    return;
  }

  const summary = document.createElement("div");
  const heading = document.createElement("h3");
  heading.textContent = tool.display_name || tool.id;
  summary.appendChild(heading);

  if (tool.description) {
    const desc = document.createElement("p");
    desc.textContent = tool.description;
    desc.style.whiteSpace = "pre-wrap";
    summary.appendChild(desc);
  }
  container.appendChild(summary);

  const form = document.createElement("form");
  form.id = "tool-form";
  form.dataset.toolId = tool.id;

  const params = Array.isArray(tool.parameters) ? tool.parameters.filter(Boolean) : [];
  const advancedParams = Array.isArray(tool.advanced_parameters)
    ? tool.advanced_parameters.filter(Boolean)
    : [];

  if (!params.length && !advancedParams.length) {
    const note = document.createElement("p");
    note.textContent = "This tool does not expose configurable parameters.";
    form.appendChild(note);
  }

  params.forEach(param => {
    form.appendChild(buildParameterField(param));
  });

  if (advancedParams.length) {
    const advanced = document.createElement("details");
    advanced.className = "advanced-parameters";

    const summary = document.createElement("summary");
    summary.textContent = "Advanced Parameters";
    advanced.appendChild(summary);

    const wrapper = document.createElement("div");
    wrapper.className = "advanced-parameters-body";
    advancedParams.forEach(param => {
      const field = buildParameterField(param);
      if (field) {
        wrapper.appendChild(field);
      }
    });

    if (wrapper.children.length) {
      advanced.appendChild(wrapper);
      form.appendChild(advanced);
    }
  }

  const submit = document.createElement("button");
  submit.type = "submit";
  submit.textContent = "Run Tool";
  form.appendChild(submit);

  form.addEventListener("submit", event => {
    event.preventDefault();
    scrollWindowToTop();
    handleToolSubmit(tool, form);
  });

  container.appendChild(form);
}

function buildInputForParameter(param, fieldKey) {
  const type = (param.type || "").toLowerCase();
  const defaultValue = param.default;
  const choiceList = Array.isArray(param.choices) ? param.choices : null;

  if (param.behavior === "dynamic") {
    const textarea = document.createElement("textarea");
    textarea.name = fieldKey;
    textarea.rows = 3;
    textarea.dataset.required = param.required ? "1" : "0";
    textarea.dataset.behavior = "dynamic";
    textarea.placeholder = "One entry per line; separate values with spaces or commas.";
    if (Array.isArray(defaultValue) && defaultValue.length) {
      const lines = defaultValue
        .map(entry => {
          if (Array.isArray(entry)) {
            return entry.join(", ");
          }
          return String(entry ?? "");
        })
        .filter(item => item.length > 0);
      textarea.value = lines.join("\n");
    }
    return textarea;
  }

  if (choiceList && choiceList.length && type !== "file" && type !== "boolean") {
    const select = document.createElement("select");
    select.name = fieldKey;
    select.dataset.required = param.required ? "1" : "0";
    select.multiple = Boolean(param.multiple);

    if (!select.multiple) {
      const placeholder = document.createElement("option");
      placeholder.value = "";
      placeholder.textContent = param.required ? "Select an option" : "Optional";
      placeholder.selected = defaultValue === undefined || defaultValue === null || defaultValue === "";
      select.appendChild(placeholder);
    }

    choiceList.forEach(choice => {
      let value = choice;
      let label = choice;
      if (choice && typeof choice === "object") {
        value = choice.value ?? choice.id ?? choice.name ?? choice.label ?? "";
        label = choice.label ?? choice.display_name ?? choice.name ?? value;
      }
      const option = document.createElement("option");
      option.value = String(value);
      option.textContent = String(label);
      select.appendChild(option);
    });

    if (select.multiple) {
      const defaults = Array.isArray(defaultValue)
        ? defaultValue.map(entry => String(entry))
        : defaultValue !== undefined && defaultValue !== null
          ? [String(defaultValue)]
          : [];
      if (defaults.length) {
        Array.from(select.options).forEach(option => {
          if (defaults.includes(option.value)) {
            option.selected = true;
          }
        });
      }
      if (choiceList.length) {
        select.size = Math.min(Math.max(choiceList.length, 2), 8);
      }
    } else if (defaultValue !== undefined && defaultValue !== null) {
      select.value = String(defaultValue);
    }

    return select;
  }

  if (type === "boolean") {
    const input = document.createElement("input");
    input.type = "checkbox";
    input.name = fieldKey;
    if (Boolean(defaultValue)) {
      input.checked = true;
    }
    input.dataset.required = param.required ? "1" : "0";
    return input;
  }

  if (
    param.multiple &&
    !choiceList &&
    type !== "file" &&
    type !== "output" &&
    type !== "boolean"
  ) {
    const textarea = document.createElement("textarea");
    textarea.name = fieldKey;
    textarea.rows = 3;
    textarea.dataset.required = param.required ? "1" : "0";
    textarea.dataset.multiline = "newline";
    if (Array.isArray(defaultValue)) {
      textarea.value = defaultValue.join("\n");
    } else if (defaultValue !== undefined && defaultValue !== null) {
      textarea.value = String(defaultValue);
    }
    return textarea;
  }

  if (type === "file") {
    const select = document.createElement("select");
    select.name = fieldKey;
    select.multiple = Boolean(param.multiple);
    select.dataset.required = param.required ? "1" : "0";
    const allowMultiple = select.multiple;
    if (!allowMultiple) {
      const placeholder = document.createElement("option");
      placeholder.value = "";
      placeholder.textContent = "Select a file";
      select.appendChild(placeholder);
    }

    const allowed = Array.isArray(param.file_types) ? param.file_types.map(v => v.toLowerCase()) : null;
    const files = allowed
      ? state.files.filter(file => file.type && allowed.includes(String(file.type).toLowerCase()))
      : state.files;

    files.forEach(file => {
      const option = document.createElement("option");
      option.value = file.file_id;
      option.textContent = `${file.original_name}${file.type ? ` (${file.type})` : ""}`;
      select.appendChild(option);
    });

    if (!files.length) {
      select.disabled = true;
      select.dataset.warning = "Upload a compatible file to use this parameter.";
    }

    if (allowMultiple && files.length) {
      select.size = Math.min(Math.max(files.length, 2), 6);
    }

    if (Array.isArray(defaultValue) && allowMultiple) {
      defaultValue.forEach(value => {
        const option = Array.from(select.options).find(opt => opt.value === value);
        if (option) {
          option.selected = true;
        }
      });
    } else if (defaultValue) {
      select.value = defaultValue;
    }

    return select;
  }

  if (type === "output") {
    const input = document.createElement("input");
    input.type = "text";
    input.name = fieldKey;
    input.placeholder = "Output filename";
    if (defaultValue) {
      input.value = defaultValue;
    }
    input.dataset.required = param.required ? "1" : "0";
    return input;
  }

  if (type === "integer") {
    const input = document.createElement("input");
    input.type = "number";
    input.name = fieldKey;
    if (defaultValue !== undefined) {
      input.value = defaultValue;
    }
    input.dataset.required = param.required ? "1" : "0";
    const stepValue = Number(param.step);
    const hasValidStep =
      param.step !== undefined &&
      param.step !== null &&
      param.step !== "" &&
      Number.isFinite(stepValue) &&
      stepValue > 0;
    input.step = hasValidStep ? String(param.step) : "1";
    input.inputMode = "numeric";
    if (param.min !== undefined && param.min !== null && param.min !== "") {
      input.min = String(param.min);
    }
    if (param.max !== undefined && param.max !== null && param.max !== "") {
      input.max = String(param.max);
    }
    return input;
  }

  if (type === "number" || type === "float") {
    const input = document.createElement("input");
    input.type = "text";
    input.name = fieldKey;
    input.autocomplete = "off";
    input.spellcheck = false;
    input.inputMode = "decimal";
    if (defaultValue !== undefined && defaultValue !== null) {
      input.value = String(defaultValue);
    }
    input.dataset.required = param.required ? "1" : "0";
    input.dataset.numericType = "float";
    const stepValue = Number(param.step);
    const hasValidStep =
      param.step !== undefined &&
      param.step !== null &&
      param.step !== "" &&
      Number.isFinite(stepValue) &&
      stepValue > 0;
    if (hasValidStep) {
      input.dataset.step = String(param.step);
    }
    if (param.min !== undefined && param.min !== null && param.min !== "") {
      input.dataset.min = String(param.min);
    }
    if (param.max !== undefined && param.max !== null && param.max !== "") {
      input.dataset.max = String(param.max);
    }
    return input;
  }

  const input = document.createElement("input");
  input.type = "text";
  input.name = fieldKey;
  if (defaultValue !== undefined) {
    input.value = defaultValue;
  }
  input.dataset.required = param.required ? "1" : "0";
  return input;
}

function buildParameterField(param) {
  const fieldKey = param.parameter || param.name;
  const displayName = param.display_name || param.name || fieldKey;

  const label = document.createElement("label");
  const name = document.createElement("span");
  name.textContent = `${displayName}${param.required ? " *" : ""}`;
  if (param.help) {
    name.title = param.help;
  }
  label.appendChild(name);

  const input = buildInputForParameter(param, fieldKey);
  if (input) {
    label.appendChild(input);
    if (input.dataset && input.dataset.warning) {
      const warn = document.createElement("span");
      warn.className = "param-help";
      warn.textContent = input.dataset.warning;
      label.appendChild(warn);
    }
  }

  if (param.help) {
    const inlineHelp = document.createElement("span");
    inlineHelp.className = "param-help";
    inlineHelp.textContent = param.help;
    label.appendChild(inlineHelp);
  }

  if (param.description) {
    const help = document.createElement("span");
    help.className = "param-help";
    help.textContent = param.description;
    label.appendChild(help);
  }

  return label;
}

function buildParameterPayload(tool, form) {
  const params = {};
  const missing = [];
  const errors = [];
  const schemaParams = [
    ...(Array.isArray(tool.parameters) ? tool.parameters : []),
    ...(Array.isArray(tool.advanced_parameters) ? tool.advanced_parameters : []),
  ];

  schemaParams.forEach(param => {
    const fieldKey = param.parameter || param.name;
    const input = form.elements.namedItem(fieldKey);
    if (!input) {
      return;
    }
    const type = (param.type || "").toLowerCase();
    let value;

    if (param.behavior === "dynamic") {
      value = parseDynamicParameter(input, param);
    } else if (type === "boolean") {
      value = input.checked;
    } else if (input instanceof HTMLSelectElement && input.multiple) {
      value = Array.from(input.selectedOptions)
        .map(option => option.value)
        .filter(Boolean);
    } else if (input instanceof HTMLInputElement || input instanceof HTMLSelectElement || input instanceof HTMLTextAreaElement) {
      value = typeof input.value === "string" ? input.value.trim() : input.value;
    } else {
      value = input.value;
    }

    const required = input.dataset.required === "1";
    if (input instanceof HTMLTextAreaElement && input.dataset.multiline === "newline") {
      value = String(value || "")
        .split(/\r?\n/)
        .map(entry => entry.trim())
        .filter(entry => entry.length > 0);
    }
    const isEmpty = value === "" || value === undefined || value === null || (Array.isArray(value) && value.length === 0);
    const displayName = param.display_name || param.name || fieldKey;
    if (required && isEmpty) {
      missing.push(displayName);
      return;
    }

    if (!required && isEmpty) {
      return;
    }

    try {
      const targetKey = param.parameter || param.name;
      params[targetKey] = coerceParameterValue(param, value);
    } catch (error) {
      errors.push(error.message);
    }
  });

  if (missing.length) {
    throw new Error(`Missing required values: ${missing.join(", ")}`);
  }
  if (errors.length) {
    throw new Error(errors.join("; "));
  }

  return params;
}

function coerceParameterValue(param, raw) {
  if (param.behavior === "dynamic") {
    return raw;
  }
  if (Array.isArray(raw)) {
    return raw.map(entry => coerceParameterValue(param, entry));
  }

  const type = (param.type || "").toLowerCase();
  const displayName = param.display_name || param.name || param.parameter || "value";
  const parseBound = bound => {
    if (bound === undefined || bound === null || bound === "") {
      return null;
    }
    const numeric = Number(bound);
    return Number.isFinite(numeric) ? numeric : null;
  };
  const minBound = parseBound(param.min);
  const maxBound = parseBound(param.max);
  const stepValue = parseBound(param.step);
  const hasStep = stepValue !== null && stepValue > 0;
  const enforceNumericConstraints = value => {
    if (minBound !== null && value < minBound) {
      throw new Error(`${displayName} must be ≥ ${minBound}`);
    }
    if (maxBound !== null && value > maxBound) {
      throw new Error(`${displayName} must be ≤ ${maxBound}`);
    }
    if (hasStep) {
      const base = minBound !== null ? minBound : 0;
      const remainder = Math.abs((value - base) / stepValue);
      const epsilon = 1e-9;
      if (Math.abs(Math.round(remainder) - remainder) > epsilon) {
        throw new Error(`${displayName} must increment by ${stepValue}`);
      }
    }
  };

  if (type === "file" || type === "output") {
    return raw;
  }

  if (type === "integer") {
    const value = Number(raw);
    if (!Number.isFinite(value)) {
      throw new Error(`${displayName} must be a number`);
    }
    if (!Number.isInteger(value)) {
      throw new Error(`${displayName} must be an integer`);
    }
    enforceNumericConstraints(value);
    return value;
  }
  if (type === "number" || type === "float") {
    const value = Number(raw);
    if (!Number.isFinite(value)) {
      throw new Error(`${displayName} must be a number`);
    }
    enforceNumericConstraints(value);
    return value;
  }
  if (type === "boolean") {
    if (typeof raw === "boolean") {
      return raw;
    }
    const lowered = String(raw).trim().toLowerCase();
    if (["true", "1", "yes", "on"].includes(lowered)) {
      return true;
    }
    if (["false", "0", "no", "off"].includes(lowered)) {
      return false;
    }
    return Boolean(raw);
  }
  return raw;
}

function parseDynamicParameter(input, param) {
  if (!(input instanceof HTMLTextAreaElement)) {
    return [];
  }
  const raw = input.value || "";
  const lines = raw
    .split(/\r?\n/)
    .map(line => line.trim())
    .filter(line => line.length > 0);

  if (!lines.length) {
    return [];
  }

  const groups = lines.map(line => line.split(/[,\s]+/).map(token => token.trim()).filter(Boolean));

  const meta = param.dynamic_action || {};
  const nargs = meta.nargs;
  const requirement = describeDynamicRequirement(nargs);
  groups.forEach((group, index) => {
    if (!validateDynamicGroup(group, nargs)) {
      const label = param.display_name || param.name || param.parameter || `parameter at line ${index + 1}`;
      throw new Error(`${label} expects ${requirement}; problem near line ${index + 1}.`);
    }
  });

  return groups;
}

function validateDynamicGroup(group, nargs) {
  if (!group.length) {
    return false;
  }
  if (nargs === undefined || nargs === null || nargs === "*") {
    return true;
  }
  if (nargs === "+") {
    return group.length >= 1;
  }
  if (nargs === "?") {
    return group.length <= 1;
  }
  const count = Number(nargs);
  if (Number.isFinite(count) && count > 0) {
    return group.length === count;
  }
  return true;
}

function describeDynamicRequirement(nargs) {
  if (nargs === undefined || nargs === null || nargs === "*") {
    return "any number of values";
  }
  if (nargs === "+") {
    return "at least one value";
  }
  if (nargs === "?") {
    return "zero or one value";
  }
  const count = Number(nargs);
  if (Number.isFinite(count) && count > 0) {
    return count === 1 ? "exactly one value" : `exactly ${count} values`;
  }
  return "valid values";
}

async function handleToolSubmit(tool, form) {
  let parameters;
  try {
    parameters = buildParameterPayload(tool, form);
  } catch (error) {
    setMessage(error.message, "error", 6000);
    return;
  }

  openPanel("jobs-panel");
  setMessage("Submitting job…", "info");
  try {
    const response = await fetch(`/api/tools/${encodeURIComponent(tool.id)}/execute`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ parameters }),
    });
    if (!response.ok) {
      const payload = await response.json().catch(() => ({}));
      const message = payload.error || `Failed to execute ${tool.id}`;
      throw new Error(message);
    }
    const result = await response.json();
    setMessage(`Job ${result.job_id} started.`, "success");
    await fetchJob(result.job_id);
  } catch (error) {
    setMessage(error.message, "error", 6000);
  }
}

async function fetchJob(jobId) {
  try {
    const response = await fetch(`/api/jobs/${encodeURIComponent(jobId)}`);
    if (!response.ok) {
      return;
    }
    const job = await response.json();
    updateJobs([job]);
  } catch (error) {
    console.error("Failed to fetch job", error);
  }
}

async function handleUploadSubmit(event) {
  event.preventDefault();
  const input = document.getElementById("file-input");
  const queuedFiles = uploadState.files.length
    ? uploadState.files
    : Array.from(input?.files || []);
  if (!queuedFiles.length) {
    setMessage("Choose one or more files first.", "error");
    return;
  }
  const formData = new FormData();
  queuedFiles.forEach(file => {
    formData.append("files", file);
  });
  setMessage("Uploading…", "info");
  try {
    const response = await fetch("/api/files/upload", {
      method: "POST",
      body: formData,
    });
    if (!response.ok) {
      const payload = await response.json().catch(() => ({}));
      const message = payload.error || "Upload failed";
      throw new Error(message);
    }
    await refreshFiles();
    clearUploadSelection();
    setMessage("Files uploaded successfully.", "success");
  } catch (error) {
    setMessage(error.message, "error", 6000);
  }
}

async function deleteFile(fileId) {
  if (!confirm("Delete this file?")) {
    return;
  }
  try {
    const response = await fetch(`/api/files/${encodeURIComponent(fileId)}`, {
      method: "DELETE",
    });
    if (!response.ok) {
      throw new Error("Failed to delete file");
    }
    setMessage("File deleted.", "success");
    await refreshFiles();
  } catch (error) {
    setMessage(error.message, "error", 6000);
  }
}

function updateJobs(jobs) {
  if (!Array.isArray(jobs) || !jobs.length) {
    renderJobs();
    return;
  }
  let latest = state.lastJobUpdate || 0;
  let shouldRefreshFiles = false;
  jobs.forEach(job => {
    state.jobs[job.job_id] = job;
    if (job.updated_at && job.updated_at > latest) {
      latest = job.updated_at;
    }
    if (
      job.status === "completed" &&
      Array.isArray(job.output_artifacts) &&
      job.output_artifacts.length
    ) {
      shouldRefreshFiles = true;
    }
  });
  if (latest > state.lastJobUpdate) {
    state.lastJobUpdate = latest;
  }
  renderJobs();
  if (shouldRefreshFiles) {
    refreshFiles();
  }
}

function renderJobs() {
  const body = document.getElementById("jobs-body");
  if (!body) {
    return;
  }
  body.innerHTML = "";
  const jobs = Object.values(state.jobs);
  if (!jobs.length) {
    const row = document.createElement("tr");
    const cell = document.createElement("td");
    cell.colSpan = 5;
    cell.textContent = "No jobs submitted yet.";
    row.appendChild(cell);
    body.appendChild(row);
    return;
  }

  jobs.sort((a, b) => (b.updated_at || b.started_at || 0) - (a.updated_at || a.started_at || 0));

  jobs.forEach(job => {
    const row = document.createElement("tr");

    const toolCell = document.createElement("td");
    toolCell.textContent = job.tool;
    row.appendChild(toolCell);

    const statusCell = document.createElement("td");
    statusCell.textContent = job.status;
    statusCell.className = `job-status ${job.status}`;
    row.appendChild(statusCell);

    const updatedCell = document.createElement("td");
    updatedCell.textContent = job.updated_at ? formatEpoch(job.updated_at) : "";
    row.appendChild(updatedCell);

    const outputsCell = document.createElement("td");
    const artifactList = Array.isArray(job.output_artifacts) ? job.output_artifacts : [];
    if (artifactList.length) {
      artifactList.forEach((artifact, index) => {
        const link = document.createElement("a");
        const artifactName = artifact.name || artifact.relative_path || artifact.file_id;
        const artifactType = artifact.type || artifact.file_type || null;
        const previewable = isPreviewableFile(artifactName, artifactType);
        if (previewable) {
          link.href = `/api/files/${encodeURIComponent(artifact.file_id)}/view`;
        } else {
          link.href = `/api/files/${encodeURIComponent(artifact.file_id)}/download`;
        }
        link.className = "download-link";
        link.textContent = artifactName;
        link.target = "_blank";
        link.rel = "noopener";
        outputsCell.appendChild(link);
        if (index < artifactList.length - 1) {
          outputsCell.appendChild(document.createElement("br"));
        }
      });
    } else if (Array.isArray(job.output_files) && job.output_files.length) {
      job.output_files.forEach((file, index) => {
        const link = document.createElement("a");
        const previewable = isPreviewableFile(file, null);
        link.href = buildJobOutputUrl(job.job_id, file, previewable ? "view" : "download");
        link.className = "download-link";
        link.textContent = file.split(/[\\/]/).pop();
        link.target = "_blank";
        link.rel = "noopener";
        outputsCell.appendChild(link);
        if (index < job.output_files.length - 1) {
          outputsCell.appendChild(document.createElement("br"));
        }
      });
    } else {
      outputsCell.textContent = job.status === "completed" ? "No outputs" : "";
    }
    row.appendChild(outputsCell);

    const configCell = document.createElement("td");
    if (job.config_file) {
      const link = document.createElement("a");
      link.href = `/api/jobs/${encodeURIComponent(job.job_id)}/config`;
      link.className = "download-link";
      link.textContent = "Download";
      link.target = "_blank";
      link.rel = "noopener";
      configCell.appendChild(link);
    }
    row.appendChild(configCell);

    const logsCell = document.createElement("td");
    if (job.log_file) {
      const link = document.createElement("a");
      link.href = `/api/jobs/${job.job_id}/logs`;
      link.textContent = "View";
      link.className = "log-link";
      link.target = "_blank";
      logsCell.appendChild(link);
    }
    if (job.status === "running" || job.status === "queued") {
      const cancelButton = document.createElement("button");
      cancelButton.type = "button";
      cancelButton.textContent = "Cancel";
      cancelButton.className = "cancel-button";
      cancelButton.setAttribute("data-cancel", job.job_id);
      logsCell.appendChild(cancelButton);
    }
    row.appendChild(logsCell);

    body.appendChild(row);
  });
}

async function cancelJob(jobId) {
  if (!jobId) {
    return;
  }
  try {
    const response = await fetch(`/api/jobs/${encodeURIComponent(jobId)}/cancel`, {
      method: "POST",
    });
    if (!response.ok) {
      const payload = await response.json().catch(() => ({}));
      const message = payload.error || "Failed to cancel job";
      throw new Error(message);
    }
    setMessage(`Job ${jobId} cancelled.`, "success");
    const job = await response.json().catch(() => null);
    if (job) {
      updateJobs([job]);
    } else {
      fetchJob(jobId);
    }
  } catch (error) {
    setMessage(error.message, "error", 6000);
  }
}

function startPolling() {
  if (pollTimer) {
    clearInterval(pollTimer);
  }
  pollTimer = setInterval(pollJobUpdates, 4000);
}

async function pollJobUpdates() {
  if (isPolling) {
    return;
  }
  isPolling = true;
  try {
    const since = state.lastJobUpdate ? `?since=${state.lastJobUpdate}` : "";
    const response = await fetch(`/api/jobs/updates${since}`);
    if (!response.ok) {
      throw new Error("Failed to fetch updates");
    }
    const jobs = await response.json();
    if (Array.isArray(jobs) && jobs.length) {
      updateJobs(jobs);
    }
  } catch (error) {
    console.error("Job polling error", error);
  } finally {
    isPolling = false;
  }
}

function setMessage(text, type = "info", timeout = 4000) {
  const box = document.getElementById("messages");
  if (!box) {
    return;
  }
  if (messageTimer) {
    clearTimeout(messageTimer);
    messageTimer = null;
  }
  box.textContent = text || "";
  box.classList.remove("info", "success", "error");
  if (text && type) {
    box.classList.add(type);
  }
  if (text && timeout) {
    messageTimer = window.setTimeout(() => {
      clearMessage();
    }, timeout);
  }
}

function clearMessage() {
  const box = document.getElementById("messages");
  if (!box) {
    return;
  }
  box.textContent = "";
  box.classList.remove("info", "success", "error");
  if (messageTimer) {
    clearTimeout(messageTimer);
    messageTimer = null;
  }
}

function formatIsoTimestamp(value) {
  if (!value) {
    return "";
  }
  const date = new Date(value);
  if (Number.isNaN(date.getTime())) {
    return value;
  }
  return date.toLocaleString();
}

function formatEpoch(seconds) {
  if (!seconds) {
    return "";
  }
  const date = new Date(seconds * 1000);
  return date.toLocaleString();
}

function formatBytes(bytes) {
  if (!bytes) {
    return "0 B";
  }
  const units = ["B", "KB", "MB", "GB", "TB"];
  const exponent = Math.min(Math.floor(Math.log(bytes) / Math.log(1024)), units.length - 1);
  const value = bytes / 1024 ** exponent;
  return `${value.toFixed(exponent === 0 ? 0 : 1)} ${units[exponent]}`;
}

function buildJobOutputUrl(jobId, relativePath, mode = "download") {
  if (!relativePath) {
    return "#";
  }
  const normalized = String(relativePath).replace(/\\+/g, "/");
  const encodedPath = normalized
    .split("/")
    .map(segment => encodeURIComponent(segment))
    .join("/");
  const base = `/api/jobs/${encodeURIComponent(jobId)}/outputs/${encodedPath}`;
  if (mode === "view") {
    return `${base}/view`;
  }
  return base;
}

function initializePanelControls() {
  const buttons = document.querySelectorAll("button[data-panel]");
  buttons.forEach(button => {
    const panelId = button.getAttribute("data-panel");
    if (!panelId) {
      return;
    }
    const panel = document.getElementById(panelId);
    if (!panel) {
      return;
    }
    button.dataset.label = button.dataset.label || button.textContent.trim();
    if (panel.classList.contains("collapsed")) {
      state.collapsedPanels.add(panelId);
    } else {
      state.collapsedPanels.delete(panelId);
    }
    updatePanelToggleButton(panelId);
  });
}

function setPanelState(panelId, open) {
  const panel = document.getElementById(panelId);
  if (!panel) {
    return;
  }
  const isOpen = !panel.classList.contains("collapsed");
  if (open === isOpen) {
    return;
  }
  if (open) {
    panel.classList.remove("collapsed");
    panel.removeAttribute("aria-hidden");
    state.collapsedPanels.delete(panelId);
  } else {
    panel.classList.add("collapsed");
    panel.setAttribute("aria-hidden", "true");
    state.collapsedPanels.add(panelId);
  }
  updatePanelToggleButton(panelId);
}

function togglePanel(panelId) {
  const shouldOpen = state.collapsedPanels.has(panelId);
  setPanelState(panelId, shouldOpen);
}

function openPanel(panelId) {
  setPanelState(panelId, true);
}

function scrollWindowToTop() {
  const prefersReducedMotion =
    window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  if (prefersReducedMotion) {
    window.scrollTo(0, 0);
  } else {
    window.scrollTo({ top: 0, behavior: "smooth" });
  }
}

function updatePanelToggleButton(panelId) {
  const button = document.querySelector(`button[data-panel="${panelId}"]`);
  if (!button) {
    return;
  }
  const isCollapsed = state.collapsedPanels.has(panelId);
  const baseLabel = button.dataset.label || button.textContent.trim() || panelId;
  button.setAttribute("aria-pressed", String(!isCollapsed));
  button.setAttribute(
    "aria-label",
    `${baseLabel} panel ${isCollapsed ? "collapsed" : "expanded"}`
  );
  button.title = `${isCollapsed ? "Show" : "Hide"} ${baseLabel}`;
  button.classList.toggle("is-collapsed", isCollapsed);
}
