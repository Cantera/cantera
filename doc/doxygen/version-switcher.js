/**
 * Version switcher for the Cantera C++ API documentation.
 *
 * Mirrors the pydata-sphinx-theme switcher used by the Sphinx docs: fetch the
 * canonical version list, populate a dropdown, and link to the same page in
 * the target version, falling back to that version's index when the page does
 * not exist there.
 *
 * If anything fails, the button keeps the static version text that Doxygen
 * substituted into it and no dropdown is attached.
 */
(function () {
    "use strict";

    // Absolute and hardcoded to /dev/, exactly as doc/sphinx/conf.py does, so
    // that released docs always show the current list of versions.
    const JSON_URL = "/dev/_static/doc-versions.json";

    /** Reduce a full version to its 'X.Y' form, as SConstruct does. */
    function shortVersion(projectNumber) {
        const match = /^\d+\.\d+/.exec(projectNumber);
        return match ? match[0] : null;
    }

    /**
     * The path of the current page relative to this version's Doxygen root.
     *
     * Doxygen is configured with CREATE_SUBDIRS, so most pages sit in a hash
     * subdirectory rather than at the root: the basename alone would not
     * resolve in another version. `$relpath^` is Doxygen's own relative path
     * from this page back to the output root, so resolving it against the
     * current URL yields that root, and the rest of the path locates the page
     * within it. The subdirectory is derived from the page name, so it agrees
     * across versions.
     */
    function currentPage(container) {
        const root = new URL(container.dataset.relpath || "./", location.href).pathname;
        if (location.pathname.startsWith(root)) {
            return location.pathname.slice(root.length) || "index.html";
        }
        return location.pathname.split("/").pop() || "index.html";
    }

    /**
     * Navigate to `href`, falling back to the target version's index page if
     * it does not exist. Older versions of Doxygen mangled page names
     * differently, so the fallback is expected to fire for 2.x and 3.0.
     */
    async function navigate(event, entry, href) {
        event.preventDefault();
        const fallback = entry.cxx_url + "index.html";
        try {
            const response = await fetch(href, {method: "HEAD"});
            location.href = response.ok ? href : fallback;
        } catch (error) {
            location.href = fallback;
        }
    }

    function populate(entries, container, button, menu) {
        const page = currentPage(container);
        const versionMatch = shortVersion(button.textContent.trim());

        entries.filter((entry) => entry && entry.cxx_url).forEach((entry) => {
            const href = entry.cxx_url + page;
            const link = document.createElement("a");
            // pydata's rule: fall back to `version` when `name` is absent.
            link.textContent = entry.name || entry.version;
            link.href = href;
            link.addEventListener("click", (event) => navigate(event, entry, href));

            if (entry.version === versionMatch) {
                link.setAttribute("aria-current", "true");
                button.textContent = link.textContent;
            }

            const item = document.createElement("li");
            item.appendChild(link);
            menu.appendChild(item);
        });

        return menu.childElementCount > 0;
    }

    function attachBehavior(container, button, menu) {
        function setExpanded(expanded) {
            button.setAttribute("aria-expanded", String(expanded));
            menu.hidden = !expanded;
        }

        button.addEventListener("click", (event) => {
            event.stopPropagation();
            setExpanded(menu.hidden);
        });

        document.addEventListener("click", (event) => {
            if (!container.contains(event.target)) {
                setExpanded(false);
            }
        });

        document.addEventListener("keydown", (event) => {
            if (event.key === "Escape") {
                setExpanded(false);
            }
        });
    }

    /**
     * Place the switcher to the right of the search box on wide screens, as it
     * sits on the rest of cantera.org, and leave it in the nav bar otherwise.
     *
     * On wide screens the switcher is moved into Doxygen's own header menu, to
     * the right of the search box. Doxygen generates that menu, so header.html
     * cannot place the element there directly -- doxygen-awesome's dark mode
     * toggle relocates itself for the same reason. The switcher and the search
     * box both float right and the first floated element lands rightmost, so
     * it is inserted *before* the search box.
     *
     * Below Doxygen's 768px breakpoint that header menu is hidden, so the
     * switcher is returned to its original place in the nav bar, where it stays
     * visible. Placement is re-evaluated when the viewport crosses the
     * breakpoint. If Doxygen's markup ever changes, the switcher stays in the
     * nav bar.
     */
    function relocate(container) {
        const wide = window.matchMedia("(min-width: 768px)");
        const navHome = container.parentNode;
        const navNext = container.nextSibling;

        function place() {
            const menu = document.getElementById("main-menu");
            const searchBox = document.getElementById("searchBoxPos2");
            if (wide.matches && menu && searchBox) {
                menu.insertBefore(container, searchBox);
            } else {
                navHome.insertBefore(container, navNext);
            }
        }

        place();
        wide.addEventListener("change", place);
    }

    async function init() {
        const container = document.getElementById("version-switcher");
        const button = document.getElementById("version-switcher-button");
        const menu = document.getElementById("version-switcher-menu");
        if (!container || !button || !menu) {
            return;
        }

        relocate(container);

        let entries;
        try {
            const response = await fetch(JSON_URL);
            if (!response.ok) {
                return;
            }
            entries = await response.json();
        } catch (error) {
            // Offline, browsing local files, or the list has moved: leave the
            // static version text in place rather than showing a broken menu.
            return;
        }

        // The list is hand-edited once per release, so treat anything that is
        // not a list of entries as a broken list rather than letting it throw:
        // an unhandled rejection here would surface as a console error.
        if (!Array.isArray(entries)) {
            return;
        }

        if (populate(entries, container, button, menu)) {
            attachBehavior(container, button, menu);
            button.classList.add("has-menu");
        }
    }

    document.addEventListener("DOMContentLoaded", init);
})();
