const state = {
  files: [],
  tools: [],
  jobs: {},
  lastJobUpdate: 0,
  activeToolId: null,
};

let messageTimer = null;
let pollTimer = null;
let isPolling = false;

document.addEventListener("DOMContentLoaded", () => {
  bindEventHandlers();
  refreshAll();
  startPolling();
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
    state.tools = Array.isArray(data) ? data : [];
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
    const downloadLink = document.createElement("a");
    downloadLink.href = `/api/files/${file.file_id}/download`;
    downloadLink.textContent = "Download";
    downloadLink.className = "download-link";
    downloadLink.target = "_blank";
    actionCell.appendChild(downloadLink);

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
      desc.textContent = tool.description;
      card.appendChild(desc);
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
    summary.appendChild(desc);
  }
  container.appendChild(summary);

  const form = document.createElement("form");
  form.id = "tool-form";
  form.dataset.toolId = tool.id;

  const params = Array.isArray(tool.parameters) ? tool.parameters : [];
  if (!params.length) {
    const note = document.createElement("p");
    note.textContent = "This tool does not expose configurable parameters.";
    form.appendChild(note);
  }

  params.forEach(param => {
    const fieldKey = param.parameter || param.name;
    const displayName = param.display_name || param.name || fieldKey;

    const label = document.createElement("label");
    const name = document.createElement("span");
    name.textContent = `${displayName}${param.required ? " *" : ""}`;
    label.appendChild(name);

    let input = buildInputForParameter(param, fieldKey);
    if (input) {
      label.appendChild(input);
      if (input.dataset && input.dataset.warning) {
        const warn = document.createElement("span");
        warn.className = "param-help";
        warn.textContent = input.dataset.warning;
        label.appendChild(warn);
      }
    }

    if (param.description) {
      const help = document.createElement("span");
      help.className = "param-help";
      help.textContent = param.description;
      label.appendChild(help);
    }

    form.appendChild(label);
  });

  const submit = document.createElement("button");
  submit.type = "submit";
  submit.textContent = "Run Tool";
  form.appendChild(submit);

  form.addEventListener("submit", event => {
    event.preventDefault();
    handleToolSubmit(tool, form);
  });

  container.appendChild(form);
}

function buildInputForParameter(param, fieldKey) {
  const type = (param.type || "").toLowerCase();
  const defaultValue = param.default;

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

  if (type === "integer" || type === "number" || type === "float") {
    const input = document.createElement("input");
    input.type = "number";
    input.name = fieldKey;
    if (defaultValue !== undefined) {
      input.value = defaultValue;
    }
    input.dataset.required = param.required ? "1" : "0";
    if (type === "integer") {
      input.step = "1";
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

function buildParameterPayload(tool, form) {
  const params = {};
  const missing = [];
  const errors = [];
  const schemaParams = Array.isArray(tool.parameters) ? tool.parameters : [];

  schemaParams.forEach(param => {
    const fieldKey = param.parameter || param.name;
    const input = form.elements.namedItem(fieldKey);
    if (!input) {
      return;
    }
    const type = (param.type || "").toLowerCase();
    let value;

    if (type === "boolean") {
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
  if (Array.isArray(raw)) {
    return raw.map(entry => coerceParameterValue(param, entry));
  }

  const type = (param.type || "").toLowerCase();

  if (type === "file" || type === "output") {
    return raw;
  }

  if (type === "integer") {
    const value = Number(raw);
    if (!Number.isFinite(value)) {
      throw new Error(`${param.display_name || param.name} must be a number`);
    }
    if (!Number.isInteger(value)) {
      throw new Error(`${param.display_name || param.name} must be an integer`);
    }
    return value;
  }
  if (type === "number" || type === "float") {
    const value = Number(raw);
    if (!Number.isFinite(value)) {
      throw new Error(`${param.display_name || param.name} must be a number`);
    }
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

async function handleToolSubmit(tool, form) {
  let parameters;
  try {
    parameters = buildParameterPayload(tool, form);
  } catch (error) {
    setMessage(error.message, "error", 6000);
    return;
  }

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
  if (!input || !input.files.length) {
    setMessage("Choose one or more files first.", "error");
    return;
  }
  const formData = new FormData();
  for (const file of input.files) {
    formData.append("files", file);
  }
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
    input.value = "";
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
        link.href = `/api/files/${encodeURIComponent(artifact.file_id)}/download`;
        link.className = "download-link";
        link.textContent = artifact.name || artifact.file_id;
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
        link.href = buildJobOutputUrl(job.job_id, file);
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

    const logsCell = document.createElement("td");
    if (job.log_file) {
      const link = document.createElement("a");
      link.href = `/api/jobs/${job.job_id}/logs`;
      link.textContent = "View";
      link.className = "log-link";
      link.target = "_blank";
      logsCell.appendChild(link);
    }
    row.appendChild(logsCell);

    body.appendChild(row);
  });
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

function buildJobOutputUrl(jobId, relativePath) {
  if (!relativePath) {
    return "#";
  }
  const normalized = String(relativePath).replace(/\\+/g, "/");
  const encodedPath = normalized
    .split("/")
    .map(segment => encodeURIComponent(segment))
    .join("/");
  return `/api/jobs/${encodeURIComponent(jobId)}/outputs/${encodedPath}`;
}
