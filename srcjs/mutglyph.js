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

let nextExportPopoverId = 0;

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
    width: 28px;
    height: 28px;
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
  .mutglyph-export-menu {
    position: relative;
    display: flex;
  }
  .mutglyph-export-popover {
    position: fixed;
    z-index: 2;
    inset: auto;
    display: grid;
    min-width: 132px;
    margin: 0;
    padding: 3px;
    border: 1px solid #bbb;
    border-radius: 4px;
    background: white;
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.18);
    color: #333;
    font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
    font-size: 12px;
  }
  .mutglyph-export-popover[hidden] {
    display: none;
  }
  .mutglyph-export-popover[popover]:not(:popover-open) {
    display: none;
  }
  @supports (position-area: bottom) {
    .mutglyph-export-popover {
      position-area: bottom span-left;
      position-try-fallbacks:
        flip-block,
        flip-inline,
        flip-block flip-inline;
      justify-self: end;
      align-self: start;
      margin-top: 4px;
    }
  }
  .mutglyph-export-option {
    padding: 6px 9px;
    border: 0;
    border-radius: 2px;
    background: transparent;
    color: #333;
    font: inherit;
    text-align: left;
    white-space: nowrap;
    cursor: pointer;
  }
  .mutglyph-export-option:hover,
  .mutglyph-export-option:focus-visible {
    background: #eee;
    outline: none;
  }
  .mutglyph-export-option:disabled {
    color: #888;
    cursor: wait;
  }
`;

function createWidgetStyle() {
  const style = document.createElement("style");
  if (typeof CSSScopeRule === "undefined") {
    // TODO: Remove this selector-prefix fallback once supported RStudio
    // Viewer versions implement prelude-less @scope.
    style.textContent = TOOLBAR_RULES.replaceAll(":scope", ".mutglyph").replace(
      /^(\s*)(\.mutglyph-)/gm,
      "$1.mutglyph $2",
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

function releasePointerFocus(button, event) {
  // Mouse clicks should not leave the subdued toolbar fully opaque. Keyboard
  // activation retains focus so the controls remain accessible without hover.
  if (event.detail > 0) queueMicrotask(() => button.blur());
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

function createImageExportOption(api, format, closeMenu) {
  const isRaster = format === "PNG";
  const title = isRaster ? "Download PNG" : "Download SVG";
  const filename = isRaster ? "mutglyph.png" : "mutglyph.svg";
  const button = document.createElement("button");
  button.type = "button";
  button.className = "mutglyph-export-option";
  button.setAttribute("role", "menuitem");
  button.textContent = title;

  button.addEventListener("click", async (event) => {
    closeMenu(event.detail > 0);
    button.disabled = true;

    try {
      const { blob, warnings = [] } = isRaster
        ? await api.imageExport.raster()
        : await api.imageExport.svg();
      for (const warning of warnings) console.warn(warning);
      downloadBlob(blob, filename);
    } catch (error) {
      console.error(`Unable to export MutGlyph as ${format}.`, error);
    } finally {
      button.disabled = false;
    }
  });

  return button;
}

function positionImageExportPopover(toggle, menu) {
  // TODO: Remove this coordinate fallback once supported RStudio Viewer
  // versions implement CSS anchor positioning.
  const gap = 4;
  const margin = 4;
  const anchor = toggle.getBoundingClientRect();
  const left = Math.min(
    anchor.right - menu.offsetWidth,
    window.innerWidth - menu.offsetWidth - margin,
  );
  let top = anchor.bottom + gap;
  if (top + menu.offsetHeight > window.innerHeight - margin) {
    top = anchor.top - menu.offsetHeight - gap;
  }
  menu.style.left = `${Math.max(margin, left)}px`;
  menu.style.top = `${Math.max(margin, top)}px`;
}

function createImageExportMenu(api) {
  const wrapper = document.createElement("div");
  wrapper.className = "mutglyph-export-menu";

  const toggle = createControlButton("image", "Download image");
  toggle.setAttribute("aria-haspopup", "menu");
  toggle.setAttribute("aria-expanded", "false");

  const menu = document.createElement("div");
  menu.className = "mutglyph-export-popover";
  menu.setAttribute("role", "menu");
  menu.setAttribute("popover", "auto");
  menu.id = `mutglyph-image-export-${++nextExportPopoverId}`;
  toggle.setAttribute("popovertarget", menu.id);

  const nativePopover = typeof menu.showPopover === "function";
  const anchorPositioning = CSS.supports?.("position-area", "bottom") ?? false;
  let dispose = () => {};
  let closeMenu;
  if (nativePopover) {
    toggle.addEventListener("click", (event) => {
      if (event.detail > 0) {
        // The popover toggle runs as the button's default action. Check its
        // final state on the next task and release focus only when it closed.
        setTimeout(() => {
          if (!menu.matches(":popover-open")) toggle.blur();
        });
      }
    });
    closeMenu = (releaseFocus = false) => {
      if (menu.matches(":popover-open")) menu.hidePopover();
      if (releaseFocus) {
        // Native popovers may restore focus to their invoker when hidden, so
        // wait until that restoration has completed before removing it.
        setTimeout(() => {
          const active = document.activeElement;
          if (active === toggle || menu.contains(active)) active.blur();
        });
      }
    };
    menu.addEventListener("toggle", (event) => {
      const open = event.newState === "open";
      toggle.setAttribute("aria-expanded", String(open));
      if (open && !anchorPositioning) positionImageExportPopover(toggle, menu);
    });
  } else {
    // Older RStudio Viewer Chromium versions lack the Popover API. Preserve
    // the same menu behavior without changing the native path above.
    toggle.removeAttribute("popovertarget");
    menu.removeAttribute("popover");
    menu.hidden = true;
    closeMenu = (releaseFocus = false) => {
      menu.hidden = true;
      toggle.setAttribute("aria-expanded", "false");
      if (releaseFocus) setTimeout(() => toggle.blur());
    };
    toggle.addEventListener("click", (event) => {
      const open = menu.hidden;
      menu.hidden = !open;
      toggle.setAttribute("aria-expanded", String(open));
      if (open) positionImageExportPopover(toggle, menu);
      else releasePointerFocus(toggle, event);
    });
    const closeOnOutsideClick = (event) => {
      if (!wrapper.contains(event.target) && !menu.contains(event.target)) {
        closeMenu();
      }
    };
    document.addEventListener("click", closeOnOutsideClick);
    dispose = () => document.removeEventListener("click", closeOnOutsideClick);
  }

  menu.append(
    createImageExportOption(api, "PNG", closeMenu),
    createImageExportOption(api, "SVG", closeMenu),
  );
  wrapper.append(toggle);

  return { element: wrapper, popover: menu, dispose };
}

function createSpecButton(spec) {
  const button = createControlButton("code", "Download GenomeSpy specification");
  button.addEventListener("click", (event) => {
    const json = JSON.stringify(spec, null, 2) + "\n";
    downloadBlob(new Blob([json], { type: "application/json" }), "mutglyph-spec.json");
    releasePointerFocus(button, event);
  });
  return button;
}

function createFullscreenButton(el) {
  const button = createControlButton("expand", "Enter fullscreen");

  if (!el.requestFullscreen) {
    button.hidden = true;
    return button;
  }

  button.addEventListener("click", async (event) => {
    releasePointerFocus(button, event);
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
  const imageExport = createImageExportMenu(api);
  controls.append(
    createFullscreenButton(el),
    imageExport.element,
    createSpecButton(spec),
  );

  return {
    element: controls,
    popover: imageExport.popover,
    dispose: imageExport.dispose,
  };
}

HTMLWidgets.widget({
  name: "mutglyph",
  type: "output",

  factory(el) {
    let api;
    let fullscreenButton;
    let disposeControls;
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

        disposeControls?.();
        disposeControls = undefined;
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
          disposeControls = controls.dispose;
          fullscreenButton = controls.element.firstElementChild;
          // Keep the popover outside the fading toolbar. A top-layer popover
          // nested under its opacity and pointer-event rules can oscillate
          // between interactive and non-interactive states during clicks.
          el.append(createWidgetStyle(), controls.element, controls.popover);
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
