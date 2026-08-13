import { embed } from "@genome-spy/core/minimal";

function fixBootstrapTooltipCollision(container) {
  // Bootstrap hides the generic `.tooltip` class used by GenomeSpy.
  // Remove this workaround when https://github.com/genome-spy/genome-spy/issues/470
  // is resolved and the bundled GenomeSpy version includes the fix.
  for (const tooltip of container.querySelectorAll(":scope > .tooltip")) {
    tooltip.style.opacity = "1";
  }
}

function createExportButton(api) {
  const button = document.createElement("button");
  button.type = "button";
  button.textContent = "Save as SVG";
  button.setAttribute("aria-label", "Save visualization as SVG");

  Object.assign(button.style, {
    position: "absolute",
    top: "8px",
    right: "8px",
    zIndex: "1",
    padding: "4px 8px",
    border: "1px solid #aaa",
    borderRadius: "3px",
    background: "rgba(255, 255, 255, 0.92)",
    color: "#333",
    font: "12px sans-serif",
    cursor: "pointer",
  });

  button.addEventListener("click", async () => {
    button.disabled = true;
    button.textContent = "Saving…";

    try {
      const { blob, warnings } = await api.imageExport.svg();
      warnings.forEach((warning) => console.warn(warning));

      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "mutglyph.svg";
      document.body.appendChild(link);
      link.click();
      link.remove();
      setTimeout(() => URL.revokeObjectURL(url), 1000);
    } catch (error) {
      console.error("Unable to export MutGlyph as SVG.", error);
    } finally {
      button.disabled = false;
      button.textContent = "Save as SVG";
    }
  });

  return button;
}

HTMLWidgets.widget({
  name: "mutglyph",
  type: "output",

  factory(el) {
    let api;
    let renderId = 0;

    el.style.position = "relative";
    el.style.overflow = "hidden";

    return {
      async renderValue(x) {
        const currentRenderId = ++renderId;

        api?.finalize();
        api = undefined;

        const container = document.createElement("div");
        container.style.width = "100%";
        container.style.height = "100%";
        el.replaceChildren(container);

        const nextApi = await embed(container, x.spec);
        fixBootstrapTooltipCollision(container);

        if (currentRenderId === renderId) {
          api = nextApi;
          el.append(createExportButton(api));
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
