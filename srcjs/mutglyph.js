import { embed } from "@genome-spy/core/minimal";

import { decodeMutGlyphTransport } from "./transport.js";

// Font Awesome Free 6.5.2 icons, https://fontawesome.com/license/free
// Only the four required SVG paths are embedded to avoid bundling its runtime,
// stylesheet, or font files. See inst/NOTICE for attribution.
const ICONS = {
  expand: {
    viewBox: "0 0 448 512",
    path: "M32 32C14.3 32 0 46.3 0 64v96c0 17.7 14.3 32 32 32s32-14.3 32-32V96h64c17.7 0 32-14.3 32-32s-14.3-32-32-32H32zM64 352c0-17.7-14.3-32-32-32s-32 14.3-32 32v96c0 17.7 14.3 32 32 32h96c17.7 0 32-14.3 32-32s-14.3-32-32-32H64V352zM320 32c-17.7 0-32 14.3-32 32s14.3 32 32 32h64v64c0 17.7 14.3 32 32 32s32-14.3 32-32V64c0-17.7-14.3-32-32-32H320zM448 352c0-17.7-14.3-32-32-32s-32 14.3-32 32v64H320c-17.7 0-32 14.3-32 32s14.3 32 32 32h96c17.7 0 32-14.3 32-32V352z",
  },
  compress: {
    viewBox: "0 0 448 512",
    path: "M160 64c0-17.7-14.3-32-32-32s-32 14.3-32 32v64H32c-17.7 0-32 14.3-32 32s14.3 32 32 32h96c17.7 0 32-14.3 32-32V64zM32 320c-17.7 0-32 14.3-32 32s14.3 32 32 32H96v64c0 17.7 14.3 32 32 32s32-14.3 32-32V352c0-17.7-14.3-32-32-32H32zM352 64c0-17.7-14.3-32-32-32s-32 14.3-32 32v96c0 17.7 14.3 32 32 32h96c17.7 0 32-14.3 32-32s-14.3-32-32-32H352V64zM320 320c-17.7 0-32 14.3-32 32v96c0 17.7 14.3 32 32 32s32-14.3 32-32V384h64c17.7 0 32-14.3 32-32s-14.3-32-32-32H320z",
  },
  image: {
    viewBox: "0 0 384 512",
    path: "M64 464c-8.8 0-16-7.2-16-16V64c0-8.8 7.2-16 16-16H224v80c0 17.7 14.3 32 32 32h80V448c0 8.8-7.2 16-16 16H64zM64 0C28.7 0 0 28.7 0 64V448c0 35.3 28.7 64 64 64H320c35.3 0 64-28.7 64-64V154.5c0-17-6.7-33.3-18.7-45.3L274.7 18.7C262.7 6.7 246.5 0 229.5 0H64zm96 256a32 32 0 1 0 -64 0 32 32 0 1 0 64 0zm69.2 46.9c-3-4.3-7.9-6.9-13.2-6.9s-10.2 2.6-13.2 6.9l-41.3 59.7-11.9-19.1c-2.9-4.7-8.1-7.5-13.6-7.5s-10.6 2.8-13.6 7.5l-40 64c-3.1 4.9-3.2 11.1-.4 16.2s8.2 8.2 14 8.2h48 32 40 72c6 0 11.4-3.3 14.2-8.6s2.4-11.6-1-16.5l-72-104z",
  },
  code: {
    viewBox: "0 0 384 512",
    path: "M64 464c-8.8 0-16-7.2-16-16V64c0-8.8 7.2-16 16-16H224v80c0 17.7 14.3 32 32 32h80V448c0 8.8-7.2 16-16 16H64zM64 0C28.7 0 0 28.7 0 64V448c0 35.3 28.7 64 64 64H320c35.3 0 64-28.7 64-64V154.5c0-17-6.7-33.3-18.7-45.3L274.7 18.7C262.7 6.7 246.5 0 229.5 0H64zm97 289c9.4-9.4 9.4-24.6 0-33.9s-24.6-9.4-33.9 0L79 303c-9.4 9.4-9.4 24.6 0 33.9l48 48c9.4 9.4 24.6 9.4 33.9 0s9.4-24.6 0-33.9l-31-31 31-31zM257 255c-9.4-9.4-24.6-9.4-33.9 0s-9.4 24.6 0 33.9l31 31-31 31c-9.4 9.4-9.4 24.6 0 33.9s24.6 9.4 33.9 0l48-48c9.4-9.4 9.4-24.6 0-33.9l-48-48z",
  },
};

const TOOLBAR_RULES = `
  .mutglyph-toolbar {
    position: absolute;
    top: 8px;
    right: 8px;
    z-index: 1;
    display: flex;
    gap: 2px;
    opacity: 0;
    pointer-events: none;
    transition: opacity 120ms ease;
  }
  :scope:hover > .mutglyph-toolbar,
  :scope:focus-within > .mutglyph-toolbar {
    opacity: 0.3;
    pointer-events: auto;
  }
  :scope > .mutglyph-toolbar:hover,
  :scope > .mutglyph-toolbar:focus-within {
    opacity: 1;
  }
  .mutglyph-control {
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 26px;
    height: 26px;
    padding: 6px;
    border: 0;
    border-radius: 3px;
    background: rgba(255, 255, 255, 0.9);
    color: #333;
    cursor: pointer;
  }
  .mutglyph-control:hover,
  .mutglyph-control:focus-visible {
    background: white;
    outline: 1px solid #999;
  }
  .mutglyph-control:disabled {
    cursor: wait;
  }
  .mutglyph-control svg {
    display: block;
    width: 100%;
    height: 100%;
    fill: currentColor;
  }
`;

function createWidgetStyle() {
  const style = document.createElement("style");
  if (typeof CSSScopeRule === "undefined") {
    // Fallback for older RStudio Viewer Chromium versions.
    style.textContent = TOOLBAR_RULES.split(":scope").join(".mutglyph").replace(
      ".mutglyph-toolbar {",
      ".mutglyph > .mutglyph-toolbar {",
    );
  } else {
    // A prelude-less scope uses the style element's parent as its root.
    style.textContent = `@scope {${TOOLBAR_RULES}}`;
  }
  return style;
}

function fixBootstrapTooltipCollision(container) {
  // Bootstrap hides the generic `.tooltip` class used by GenomeSpy.
  // Remove this workaround when https://github.com/genome-spy/genome-spy/issues/470
  // is resolved and the bundled GenomeSpy version includes the fix.
  for (const tooltip of container.querySelectorAll(":scope > .tooltip")) {
    tooltip.style.opacity = "1";
  }
}

function createIcon(name) {
  const icon = ICONS[name];
  const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
  svg.setAttribute("viewBox", icon.viewBox);
  svg.setAttribute("aria-hidden", "true");
  path.setAttribute("d", icon.path);
  svg.append(path);
  return svg;
}

function setButton(button, icon, title) {
  button.replaceChildren(createIcon(icon));
  button.title = title;
  button.setAttribute("aria-label", title);
}

function createControlButton(icon, title) {
  const button = document.createElement("button");
  button.type = "button";
  button.className = "mutglyph-control";
  setButton(button, icon, title);
  return button;
}

function downloadBlob(blob, filename) {
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function createExportButton(api) {
  const button = createControlButton("image", "Download SVG");

  button.addEventListener("click", async () => {
    button.disabled = true;

    try {
      const { blob, warnings } = await api.imageExport.svg();
      warnings.forEach((warning) => console.warn(warning));
      downloadBlob(blob, "mutglyph.svg");
    } catch (error) {
      console.error("Unable to export MutGlyph as SVG.", error);
    } finally {
      button.disabled = false;
    }
  });

  return button;
}

function createSpecButton(spec) {
  const button = createControlButton("code", "Download GenomeSpy specification");
  button.addEventListener("click", () => {
    const json = JSON.stringify(spec, null, 2) + "\n";
    downloadBlob(new Blob([json], { type: "application/json" }), "mutglyph-spec.json");
  });
  return button;
}

function createFullscreenButton(el) {
  const button = createControlButton("expand", "Enter fullscreen");

  if (!el.requestFullscreen) {
    button.hidden = true;
    return button;
  }

  button.addEventListener("click", async () => {
    try {
      if (document.fullscreenElement === el) {
        await document.exitFullscreen();
      } else {
        await el.requestFullscreen();
      }
    } catch (error) {
      console.error("Unable to change MutGlyph fullscreen mode.", error);
    }
  });

  return button;
}

function createControls(api, el, spec) {
  const controls = document.createElement("div");
  controls.className = "mutglyph-toolbar";
  controls.append(
    createFullscreenButton(el),
    createExportButton(api),
    createSpecButton(spec),
  );
  return controls;
}

HTMLWidgets.widget({
  name: "mutglyph",
  type: "output",

  factory(el) {
    let api;
    let fullscreenButton;
    let renderId = 0;

    el.style.position = "relative";
    el.style.overflow = "hidden";
    el.addEventListener("fullscreenchange", () => {
      if (!fullscreenButton) return;
      const active = document.fullscreenElement === el;
      setButton(
        fullscreenButton,
        active ? "compress" : "expand",
        active ? "Exit fullscreen" : "Enter fullscreen",
      );
    });

    return {
      async renderValue(x) {
        const currentRenderId = ++renderId;
        const spec = decodeMutGlyphTransport(x.spec);

        api?.finalize();
        api = undefined;

        const container = document.createElement("div");
        container.style.width = "100%";
        container.style.height = "100%";
        el.replaceChildren(container);

        const nextApi = await embed(container, spec);
        fixBootstrapTooltipCollision(container);

        if (currentRenderId === renderId) {
          api = nextApi;
          const controls = createControls(api, el, spec);
          fullscreenButton = controls.firstElementChild;
          el.append(createWidgetStyle(), controls);
        } else {
          nextApi.finalize();
          container.remove();
        }
      },

      resize() {
        // GenomeSpy observes its container with ResizeObserver.
      },
    };
  },
});
